#!/usr/bin/env bash
# Reproducible HaplotypeCaller Peak-RSS profile: gatk-rs (Rust) vs pinned Java GATK 4.4.
#
# Runs two labeled profiles by default:
#   smoke      — checked-in 32 bp fixture (runtime-overhead reference only)
#   realistic  — multi-Mb GIAB-dense NA12878 window (only this may back public
#                "X% less memory" claims, and only when measured on the dedicated
#                gatk-rs-benchmark host — see docs/ci/PERF_BENCHMARK_HOST.md)
#
# Output:
#   docs/perf/HC_MEMORY_PROFILE.md          — human-readable report
#   docs/perf/hc_memory_profile_latest.json — machine-readable summary
#   docs/perf/runs/<stamp>/                — raw time logs + VCFs
#
# Requirements:
#   - cargo, samtools
#   - Java GATK via GATK_JAR, `gatk` on PATH, or Docker image from docs/GATK_PINNED.env
#   - macOS: /usr/bin/time -l ; Linux: /usr/bin/time -v (or gtime)
#   - realistic profile: network + disk to stage GIAB BAM slice (unless pre-staged)
#
# Usage:
#   ./scripts/perf/run_hc_memory_profile.sh
#   HC_MEM_PROFILES=smoke ./scripts/perf/run_hc_memory_profile.sh
#   HC_MEM_PROFILES=realistic JAVA_XMX=4g ./scripts/perf/run_hc_memory_profile.sh
#
# IMPORTANT: Do not invent README "X% less memory" claims from the smoke profile.
# Public memory claims require the realistic profile on the dedicated bench host.
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"
# shellcheck source=../parity/lib_pinned_gatk.sh
source "${repo_root}/scripts/parity/lib_pinned_gatk.sh"
# shellcheck source=../parity/giab/lib_giab.sh
source "${repo_root}/scripts/parity/giab/lib_giab.sh"

stamp="$(date -u +%Y%m%dT%H%M%SZ)"
out_root="${repo_root}/docs/perf"
run_dir="${out_root}/runs/${stamp}"
mkdir -p "${run_dir}"

# Comma-separated: smoke,realistic
HC_MEM_PROFILES="${HC_MEM_PROFILES:-smoke,realistic}"

# --- Smoke fixture (tiny, checked-in) ---
SMOKE_REF="${HC_MEM_REF:-${repo_root}/parity/fixtures/reference.fa}"
SMOKE_BAM="${HC_MEM_BAM:-${repo_root}/parity/fixtures/sample.bam}"
SMOKE_INTERVAL="${HC_MEM_INTERVAL:-chr1:1-32}"

# --- Realistic fixture (GIAB-dense multi-Mb window; separate OUT_DIR from 50kb fair-HC slice) ---
# Default: 2 Mb on chr20 starting at the known-dense 10 Mb locus used elsewhere in the repo.
REAL_INTERVAL="${HC_MEM_REALISTIC_INTERVAL:-20:10000000-12000000}"
REAL_OUT_DIR="${HC_MEM_REALISTIC_OUT_DIR:-${repo_root}/parity/realworld/na12878_giab_window_mem_2mb_b37}"
REAL_REF="${HC_MEM_REALISTIC_REF:-${repo_root}/parity/realworld/assets/hs37d5.simple.fa}"
REAL_BAM="${HC_MEM_REALISTIC_BAM:-${REAL_OUT_DIR}/NA12878_giab_window.b37.bam}"
HC_MEM_STAGE_GIAB="${HC_MEM_STAGE_GIAB:-1}"

# Realistic pipeline heap for GATK HC shards (Broad WDLs commonly use 4g–6g).
JAVA_XMX="${JAVA_XMX:-4g}"
JAVA_XMS="${JAVA_XMS:-1g}"
JAVA_OPTS="${JAVA_OPTS:--Xms${JAVA_XMS} -Xmx${JAVA_XMX}}"
# Match fair-HC default: single-threaded Peak-RSS (multi-thread inflates RSS / can OOM).
HC_MEM_THREADS="${HC_MEM_THREADS:-1}"
export RAYON_NUM_THREADS="${RAYON_NUM_THREADS:-${HC_MEM_THREADS}}"

export CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-1}"
target_dir="${CARGO_TARGET_DIR:-${repo_root}/target}"
echo "[hc-mem] building release gatk-cli (CARGO_TARGET_DIR=${target_dir})…"
echo "[hc-mem] HC_MEM_THREADS=${HC_MEM_THREADS}"
(
  cd "${repo_root}"
  cargo build -p gatk-cli --release --locked
)
RUST_BIN="${target_dir}/release/gatk-rs"
if [[ ! -x "${RUST_BIN}" ]]; then
  echo "missing ${RUST_BIN}" >&2
  exit 1
fi

backend="$(giab_time_backend)"
echo "[hc-mem] time backend=${backend} profiles=${HC_MEM_PROFILES}"

stage_realistic_inputs() {
  if [[ ! -f "${REAL_REF}" ]]; then
    echo "[hc-mem] staging hs37d5.simple.fa…"
    "${repo_root}/scripts/parity/realworld/03_stage_reference_and_truth.sh"
  fi
  if [[ -f "${REAL_BAM}" && ( -f "${REAL_BAM}.bai" || -f "${REAL_BAM%.*}.bai" ) ]]; then
    local staged_iv=""
    if [[ -f "${REAL_OUT_DIR}/stage_meta.json" ]]; then
      staged_iv="$(python3 -c "import json; print(json.load(open('${REAL_OUT_DIR}/stage_meta.json')).get('interval',''))" 2>/dev/null || true)"
    fi
    if [[ "${staged_iv}" == "${REAL_INTERVAL}" ]]; then
      echo "[hc-mem] reuse staged realistic BAM for ${REAL_INTERVAL}"
      return 0
    fi
  fi
  if [[ "${HC_MEM_STAGE_GIAB}" != "1" ]]; then
    echo "[hc-mem] ERROR: missing realistic BAM ${REAL_BAM} (set HC_MEM_STAGE_GIAB=1)" >&2
    exit 1
  fi
  echo "[hc-mem] staging GIAB-dense BAM slice ${REAL_INTERVAL} → ${REAL_OUT_DIR}"
  J6_DENSE_INTERVAL="${REAL_INTERVAL}" \
    J6_DENSE_OUT_DIR="${REAL_OUT_DIR}" \
    "${repo_root}/scripts/parity/realworld/04_stage_na12878_giab_window_bam.sh"
}

# Run one HC pair (Rust + Java); writes time logs + summary JSON fragment to stdout path.
# Args: profile_id label ref bam interval note marketing_ok(0|1)
run_profile_pair() {
  local profile_id="$1" label="$2" ref="$3" bam="$4" interval="$5" note="$6" marketing_ok="$7"
  local prof_dir="${run_dir}/${profile_id}"
  mkdir -p "${prof_dir}"

  if [[ ! -f "${bam}.bai" && ! -f "${bam%.*}.bai" ]]; then
    samtools index "${bam}"
  fi

  local rust_vcf="${prof_dir}/rust.hc.vcf"
  local java_vcf="${prof_dir}/java.hc.vcf"
  local rust_time="${prof_dir}/rust.time.txt"
  local java_time="${prof_dir}/java.time.txt"
  local rust_log="${prof_dir}/rust.stdout.txt"
  local java_log="${prof_dir}/java.stdout.txt"
  local java_cmd_file="${prof_dir}/java.cmdline.txt"
  local summary_json="${prof_dir}/summary.json"

  echo "[hc-mem] === profile=${profile_id} (${label}) interval=${interval} threads=${HC_MEM_THREADS} ==="

  echo "[hc-mem] Rust HaplotypeCaller…"
  set +e
  giab_run_timed "${backend}" "${rust_time}" \
    "${RUST_BIN}" HaplotypeCaller \
      -R "${ref}" \
      -I "${bam}" \
      -O "${rust_vcf}" \
      -L "${interval}" \
      --threads "${HC_MEM_THREADS}" \
    >"${rust_log}" 2>&1
  local rust_exit=$?
  set -e
  if [[ -f "${rust_time}.stdout" ]]; then
    cat "${rust_time}.stdout" >>"${rust_log}" || true
  fi
  if [[ "${rust_exit}" -ne 0 ]]; then
    echo "Rust HaplotypeCaller failed (profile=${profile_id} exit=${rust_exit}). See ${rust_log}" >&2
    # macOS /usr/bin/time often prints the real failure on its stderr log.
    tail -40 "${rust_log}" >&2 || true
    tail -20 "${rust_time}" >&2 || true
    return "${rust_exit}"
  fi

  echo "[hc-mem] Java GATK HaplotypeCaller (${GATK_DOCKER_IMAGE:-local})…"
  run_java_timed() {
    if [[ -n "${GATK_JAR:-}" && -f "${GATK_JAR}" ]]; then
      echo "java ${JAVA_OPTS} -jar ${GATK_JAR} HaplotypeCaller -R ${ref} -I ${bam} -O ${java_vcf} -L ${interval} --native-pair-hmm-threads ${HC_MEM_THREADS}" \
        | tee "${java_cmd_file}"
      giab_run_timed "${backend}" "${java_time}" \
        java ${JAVA_OPTS} -jar "${GATK_JAR}" HaplotypeCaller \
          -R "${ref}" -I "${bam}" -O "${java_vcf}" -L "${interval}" \
          --native-pair-hmm-threads "${HC_MEM_THREADS}" \
        >"${java_log}" 2>&1
      return $?
    fi
    if command -v gatk >/dev/null 2>&1; then
      echo "gatk --java-options '${JAVA_OPTS}' HaplotypeCaller -R ${ref} -I ${bam} -O ${java_vcf} -L ${interval} --native-pair-hmm-threads ${HC_MEM_THREADS}" \
        | tee "${java_cmd_file}"
      giab_run_timed "${backend}" "${java_time}" \
        gatk --java-options "${JAVA_OPTS}" HaplotypeCaller \
          -R "${ref}" -I "${bam}" -O "${java_vcf}" -L "${interval}" \
          --native-pair-hmm-threads "${HC_MEM_THREADS}" \
        >"${java_log}" 2>&1
      return $?
    fi
    if [[ -z "${GATK_DOCKER_IMAGE:-}" ]]; then
      echo "No Java GATK available" >&2
      return 127
    fi
    local host_root="${PARITY_HOST_REPO_ROOT:-${repo_root}}"
    local plat_args=()
    [[ -n "${GATK_DOCKER_PLATFORM:-}" ]] && plat_args+=(--platform "${GATK_DOCKER_PLATFORM}")
    {
      echo "docker run --rm ${plat_args[*]:-} -v ${host_root}:${host_root} -w ${host_root} ${GATK_DOCKER_IMAGE} \\"
      echo "  # Peak RSS via /proc VmHWM sampler (image has no GNU time)"
      echo "  gatk --java-options '${JAVA_OPTS}' HaplotypeCaller \\"
      echo "  -R ${ref} -I ${bam} -O ${java_vcf} -L ${interval} --native-pair-hmm-threads ${HC_MEM_THREADS}"
    } | tee "${java_cmd_file}"

    local inner="${prof_dir}/java_inner.sh"
    cat >"${inner}" <<EOF
#!/bin/bash
set -euo pipefail
gatk --java-options '${JAVA_OPTS}' HaplotypeCaller \\
  -R '${ref}' -I '${bam}' -O '${java_vcf}' -L '${interval}' \\
  --native-pair-hmm-threads '${HC_MEM_THREADS}' \\
  >'${java_log}' 2>&1 &
gp=\$!
start_ns=\$(date +%s%N)
peak_kb=0
while kill -0 "\$gp" 2>/dev/null; do
  for status in /proc/[0-9]*/status; do
    [[ -r "\$status" ]] || continue
    name=\$(awk '/^Name:/ {print \$2}' "\$status" 2>/dev/null || true)
    case "\$name" in
      java|gatk) ;;
      *) continue ;;
    esac
    hwm=\$(awk '/^VmHWM:/ {print \$2}' "\$status" 2>/dev/null || true)
    [[ -n "\$hwm" ]] || continue
    if (( hwm > peak_kb )); then peak_kb=\$hwm; fi
  done
  sleep 0.05
done
wait "\$gp"
ec=\$?
end_ns=\$(date +%s%N)
wall=\$(awk -v s="\$start_ns" -v e="\$end_ns" 'BEGIN { printf "%.3f", (e - s) / 1e9 }')
{
  echo "Elapsed (wall clock) time (h:mm:ss or m:ss): 0:\${wall}"
  echo "Maximum resident set size (kbytes): \${peak_kb}"
  echo "Command being timed: gatk HaplotypeCaller (Linux /proc VmHWM sampler; image lacks GNU time)"
} >'${java_time}'
exit "\$ec"
EOF
    chmod +x "${inner}"
    docker run --rm \
      "${plat_args[@]}" \
      -v "${host_root}:${host_root}" \
      -w "${host_root}" \
      "${GATK_DOCKER_IMAGE}" \
      /bin/bash "${inner}"
  }

  set +e
  run_java_timed
  local java_exit=$?
  set -e
  if [[ "${java_exit}" -ne 0 ]]; then
    echo "Java HaplotypeCaller failed (profile=${profile_id} exit=${java_exit}). See ${java_log}" >&2
    tail -60 "${java_log}" >&2 || true
    return "${java_exit}"
  fi

  local rust_json java_json
  rust_json="$(giab_parse_time_log "${rust_time}")"
  java_json="$(giab_parse_time_log "${java_time}")"

  python3 - "${summary_json}" "${profile_id}" "${label}" \
    "${ref}" "${bam}" "${interval}" "${note}" "${marketing_ok}" \
    "${JAVA_OPTS}" "${GATK_DOCKER_IMAGE:-}" "${GATK_PINNED_SHA:-}" \
    "${rust_json}" "${java_json}" "${RUST_BIN}" <<'PY'
import json, sys, pathlib

(
    out_path, profile_id, label,
    ref, bam, interval, note, marketing_ok,
    java_opts, docker_image, gatk_sha,
    rust_json, java_json, rust_bin,
) = sys.argv[1:]

rust = json.loads(rust_json)
java = json.loads(java_json)
rk, jk = rust.get("max_rss_kb"), java.get("max_rss_kb")
summary = {
    "profile_id": profile_id,
    "label": label,
    "marketing_claim_ok": marketing_ok == "1",
    "fixture": {
        "reference": ref,
        "bam": bam,
        "interval": interval,
        "note": note,
    },
    "rust": {
        "binary": rust_bin,
        "wall_sec": rust.get("wall_sec"),
        "max_rss_kb": rk,
        "max_rss_mib": None if rk is None else rk / 1024.0,
        "threads": int(__import__("os").environ.get("HC_MEM_THREADS", "1")),
        "cmdline": (
            f"{rust_bin} HaplotypeCaller -R {ref} -I {bam} "
            f"-O <run>/{profile_id}/rust.hc.vcf -L {interval} "
            f"--threads {__import__('os').environ.get('HC_MEM_THREADS', '1')}"
        ),
    },
    "java": {
        "gatk_version_pin": "4.4.0.0",
        "gatk_sha": gatk_sha,
        "docker_image": docker_image or None,
        "java_options": java_opts,
        "native_pair_hmm_threads": int(__import__("os").environ.get("HC_MEM_THREADS", "1")),
        "wall_sec": java.get("wall_sec"),
        "max_rss_kb": jk,
        "max_rss_mib": None if jk is None else jk / 1024.0,
        "note": "Peak RSS = max /proc VmHWM (kB) of java|gatk inside Docker when the image lacks GNU time; else GNU time -v on the JVM process.",
    },
}
if rk and jk and jk > 0:
    summary["ratio"] = {
        "java_over_rust_peak_rss": jk / rk,
        "rust_fraction_of_java_peak_rss": rk / jk,
        "delta_kb": jk - rk,
    }
pathlib.Path(out_path).write_text(json.dumps(summary, indent=2) + "\n")
print(f"[hc-mem] wrote {out_path}  Rust={rk} KiB  Java={jk} KiB")
PY
}

write_failed_profile_stub() {
  local profile_id="$1" label="$2" ref="$3" bam="$4" interval="$5" note="$6" marketing_ok="$7" err="$8"
  local prof_dir="${run_dir}/${profile_id}"
  mkdir -p "${prof_dir}"
  python3 - "${prof_dir}/summary.json" "${profile_id}" "${label}" \
    "${ref}" "${bam}" "${interval}" "${note}" "${marketing_ok}" "${err}" <<'PY'
import json, pathlib, sys
out, pid, label, ref, bam, interval, note, m_ok, err = sys.argv[1:]
pathlib.Path(out).write_text(json.dumps({
    "profile_id": pid,
    "label": label,
    "marketing_claim_ok": m_ok == "1",
    "status": "failed",
    "error": err,
    "fixture": {"reference": ref, "bam": bam, "interval": interval, "note": note},
    "rust": {"wall_sec": None, "max_rss_kb": None, "max_rss_mib": None},
    "java": {"wall_sec": None, "max_rss_kb": None, "max_rss_mib": None,
             "gatk_version_pin": "4.4.0.0"},
}, indent=2) + "\n")
PY
}

# --- Select and run profiles ---
IFS=',' read -r -a profile_list <<<"${HC_MEM_PROFILES}"
declare -a ran_profiles=()
profile_failures=0
for raw in "${profile_list[@]}"; do
  p="$(echo "${raw}" | tr -d '[:space:]')"
  [[ -n "${p}" ]] || continue
  case "${p}" in
    smoke)
      if run_profile_pair \
        "smoke" \
        "Trivial smoke (reproducibility only)" \
        "${SMOKE_REF}" "${SMOKE_BAM}" "${SMOKE_INTERVAL}" \
        "Checked-in p4 smoke fixture (chr1 32bp). Peak RSS is dominated by process/runtime overhead; never use for public memory claims." \
        "0"
      then
        ran_profiles+=("smoke")
      else
        echo "[hc-mem] ERROR: smoke profile is required and failed" >&2
        exit 1
      fi
      ;;
    realistic)
      stage_realistic_inputs
      if run_profile_pair \
        "realistic" \
        "Realistic GIAB-dense multi-Mb window" \
        "${REAL_REF}" "${REAL_BAM}" "${REAL_INTERVAL}" \
        "NA12878 GIAB-dense window (default 20:10000000-12000000, 2 Mb) staged from NIST 30× BAM. Suitable as the basis for public memory claims only when measured on the dedicated gatk-rs-benchmark host with HOST_SPECS.md populated." \
        "1"
      then
        ran_profiles+=("realistic")
      else
        echo "[hc-mem] WARN: realistic profile failed — writing stub (re-run on dedicated gatk-rs-benchmark host)" >&2
        write_failed_profile_stub \
          "realistic" \
          "Realistic GIAB-dense multi-Mb window" \
          "${REAL_REF}" "${REAL_BAM}" "${REAL_INTERVAL}" \
          "NA12878 GIAB-dense window (default 20:10000000-12000000, 2 Mb). Measurement failed on this host (often OOM on laptops); re-run on dedicated gatk-rs-benchmark." \
          "1" \
          "engine exited non-zero (see docs/perf/runs/${stamp}/realistic/)"
        ran_profiles+=("realistic")
        profile_failures=$((profile_failures + 1))
      fi
      ;;
    *)
      echo "[hc-mem] unknown profile '${p}' (expected smoke|realistic)" >&2
      exit 2
      ;;
  esac
done

if [[ "${#ran_profiles[@]}" -eq 0 ]]; then
  echo "[hc-mem] ERROR: no profiles ran (HC_MEM_PROFILES=${HC_MEM_PROFILES})" >&2
  exit 2
fi

rustc_version="$(rustc --version 2>/dev/null || echo unknown)"
cargo_version="$(cargo --version 2>/dev/null || echo unknown)"
uname_s="$(uname -srm)"
rust_sha="$(git -C "${repo_root}" rev-parse --short HEAD 2>/dev/null || echo unknown)"
# Dedicated-host gate: publishable only when HOST_SPECS was captured for this runner label.
dedicated_ok=0
if [[ -f "${out_root}/HOST_SPECS.md" ]] \
  && grep -q 'gatk-rs-benchmark' "${out_root}/HOST_SPECS.md" 2>/dev/null \
  && ! grep -q 'not yet captured' "${out_root}/HOST_SPECS.md" 2>/dev/null; then
  dedicated_ok=1
fi

python3 - "${run_dir}" "${out_root}" "${stamp}" \
  "${JAVA_OPTS}" "${GATK_DOCKER_IMAGE:-}" "${GATK_PINNED_SHA:-}" \
  "${rustc_version}" "${cargo_version}" "${uname_s}" "${rust_sha}" \
  "${RUST_BIN}" "${dedicated_ok}" "${ran_profiles[*]}" <<'PY'
import json, pathlib, sys

(
    run_dir, out_root, stamp,
    java_opts, docker_image, gatk_sha,
    rustc_version, cargo_version, uname_s, rust_sha,
    rust_bin, dedicated_ok, ran_profiles_s,
) = sys.argv[1:]

ran = ran_profiles_s.split()
profiles = []
for pid in ran:
    p = pathlib.Path(run_dir) / pid / "summary.json"
    profiles.append(json.loads(p.read_text()))

def fmt_mib(kb):
    if kb is None:
        return "n/a"
    return f"{kb / 1024.0:.2f} MiB ({kb} KiB)"

dedicated = dedicated_ok == "1"
by_id = {p["profile_id"]: p for p in profiles}

# Redact machine-local absolute paths before committing docs.
repo = str(pathlib.Path(out_root).resolve().parent.parent)
def redact(obj):
    if isinstance(obj, str):
        return obj.replace(repo, "<repo>").replace(rust_bin, "<release-bin>/gatk-rs")
    if isinstance(obj, list):
        return [redact(x) for x in obj]
    if isinstance(obj, dict):
        return {k: redact(v) for k, v in obj.items()}
    return obj
profiles = redact(profiles)
by_id = {p["profile_id"]: p for p in profiles}

summary = {
    "stamp_utc": stamp,
    "host": uname_s,
    "git_sha": rust_sha,
    "rustc": rustc_version,
    "cargo": cargo_version,
    "binary": "<release-bin>/gatk-rs",
    "build": "cargo build -p gatk-cli --release --locked",
    "java_options": java_opts,
    "gatk_docker_image": docker_image or None,
    "gatk_sha": gatk_sha,
    "dedicated_benchmark_host": dedicated,
    "public_memory_claim_allowed": bool(
        dedicated
        and "realistic" in by_id
        and by_id["realistic"].get("status") != "failed"
        and by_id["realistic"].get("rust", {}).get("max_rss_kb")
        and by_id["realistic"].get("java", {}).get("max_rss_kb")
    ),
    "profiles": profiles,
}

# Backward-compatible top-level smoke fields (older consumers).
if "smoke" in by_id:
    sm = by_id["smoke"]
    summary["fixture"] = sm["fixture"]
    summary["rust"] = {**sm["rust"], "git_sha": rust_sha, "build": summary["build"],
                       "rustc": rustc_version, "cargo": cargo_version}
    summary["java"] = sm["java"]
    if "ratio" in sm:
        summary["ratio"] = sm["ratio"]

pathlib.Path(out_root).mkdir(parents=True, exist_ok=True)
json_path = pathlib.Path(out_root) / "hc_memory_profile_latest.json"
json_path.write_text(json.dumps(summary, indent=2) + "\n")
(pathlib.Path(run_dir) / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

def profile_table(p):
    if p.get("status") == "failed":
        return (
            f"**Status:** measurement **failed** on this host "
            f"(`{p.get('error', 'unknown')}`).\n\n"
            "Re-run on the dedicated `gatk-rs-benchmark` host "
            "([PERF_BENCHMARK_HOST.md](../ci/PERF_BENCHMARK_HOST.md)).\n"
        )
    rk = p["rust"].get("max_rss_kb")
    jk = p["java"].get("max_rss_kb")
    rows = f"""| Engine | Peak RSS | Wall time |
|--------|----------|-----------|
| **gatk-rs** (Rust release) | **{fmt_mib(rk)}** | {p['rust'].get('wall_sec') if p['rust'].get('wall_sec') is not None else 'n/a'} s |
| **Java GATK 4.4.0.0** | **{fmt_mib(jk)}** | {p['java'].get('wall_sec') if p['java'].get('wall_sec') is not None else 'n/a'} s |
"""
    if "ratio" in p:
        r = p["ratio"]
        rows += f"""
| Java / Rust Peak-RSS | {r['java_over_rust_peak_rss']:.2f}× |
| Rust as fraction of Java Peak-RSS | {100.0 * r['rust_fraction_of_java_peak_rss']:.1f}% |
| Absolute delta (Java − Rust) | {r['delta_kb']/1024.0:.2f} MiB |
"""
    return rows

sections = []
if "smoke" in by_id:
    p = by_id["smoke"]
    sections.append(f"""## A. Trivial smoke — reproducibility reference only

**Label:** {p['label']}  
**Interval:** `{p['fixture']['interval']}`  
**BAM / ref:** checked-in `parity/fixtures/`  

> **Not for marketing.** Peak-RSS here is dominated by JVM/runtime fixed cost
> on a 32 bp window. Do **not** derive “X% less memory” from this table.

{profile_table(p)}
""")

if "realistic" in by_id:
    p = by_id["realistic"]
    claim = (
        "This realistic profile was measured on a host that reports "
        "`gatk-rs-benchmark` in [`HOST_SPECS.md`](HOST_SPECS.md) — it **may** "
        "back a public memory claim (cite this report + HOST_SPECS)."
        if dedicated
        else "This realistic profile was **not** measured on the dedicated "
        "`gatk-rs-benchmark` host (see [`HOST_SPECS.md`](HOST_SPECS.md) / "
        "[`docs/ci/PERF_BENCHMARK_HOST.md`](../ci/PERF_BENCHMARK_HOST.md)). "
        "Numbers below are engineering evidence only — **do not** use them for "
        "a public “X% less memory” claim until re-run on that host."
    )
    sections.append(f"""## B. Realistic GIAB-dense window — public-claim basis

**Label:** {p['label']}  
**Interval:** `{p['fixture']['interval']}` (multi-Mb; default 2 Mb on chr20 dense locus)  
**BAM:** staged NA12878 NIST 30× slice (`parity/realworld/na12878_giab_window_mem_2mb_b37/`)  

> {claim}

{profile_table(p)}
""")

claim_line = (
    "**Public memory claim status:** allowed from profile **B** on this host."
    if summary["public_memory_claim_allowed"]
    else "**Public memory claim status:** **not allowed** from this run "
    "(need realistic profile on dedicated `gatk-rs-benchmark` host)."
)

md = f"""# HaplotypeCaller memory profile (reproducible)

**Generated (UTC):** `{stamp}`  
**Host:** `{uname_s}`  
**Git:** `{rust_sha}`  
**Runner script:** [`scripts/perf/run_hc_memory_profile.sh`](../../scripts/perf/run_hc_memory_profile.sh)  
**Raw run directory:** `docs/perf/runs/{stamp}/`

{claim_line}

Profiles measured: `{', '.join(ran)}`.

{''.join(sections)}
## Exact commands

### Rust (build once)

```bash
cargo build -p gatk-cli --release --locked
# rustc: {rustc_version}
# cargo: {cargo_version}
# git: {rust_sha}
```

Per-profile command lines are in each `docs/perf/runs/{stamp}/<profile>/` summary
and the JSON under `rust.cmdline` / `java.cmdline.txt`.

### Java GATK 4.4

- Pin: `GATK_PINNED_SHA={gatk_sha}` (`docs/GATK_PINNED.env`)
- Image: `{docker_image or 'local gatk / GATK_JAR'}`
- JVM options (pipeline-realistic): `{java_opts}`

```bash
./scripts/perf/run_hc_memory_profile.sh
# smoke only:
#   HC_MEM_PROFILES=smoke ./scripts/perf/run_hc_memory_profile.sh
# realistic only (stages 2 Mb GIAB window if needed):
#   HC_MEM_PROFILES=realistic ./scripts/perf/run_hc_memory_profile.sh
```

When Docker is used, Peak-RSS is sampled from `/proc/*/status` **VmHWM**
for `java`/`gatk` **inside** the Linux container (the Broad 4.4 image has no
GNU `/usr/bin/time`). Host `time docker …` is never used for RSS.

## Re-run

```bash
./scripts/perf/run_hc_memory_profile.sh
# overrides:
#   HC_MEM_PROFILES=smoke,realistic
#   HC_MEM_REALISTIC_INTERVAL=20:10000000-12000000
#   JAVA_XMX=4g JAVA_XMS=1g
# Dedicated host (published claims):
#   see docs/ci/PERF_BENCHMARK_HOST.md + Actions workflow benchmark.yml
```
"""
md_path = pathlib.Path(out_root) / "HC_MEMORY_PROFILE.md"
md_path.write_text(md)
print(f"Wrote {md_path}")
print(f"Wrote {json_path}")
print(f"dedicated_benchmark_host={dedicated} public_memory_claim_allowed={summary['public_memory_claim_allowed']}")
for p in profiles:
    print(
        f"  {p['profile_id']}: rust_kib={p['rust'].get('max_rss_kb')} "
        f"java_kib={p['java'].get('max_rss_kb')}"
    )
PY

echo "[hc-mem] done → docs/perf/HC_MEMORY_PROFILE.md"
if [[ "${profile_failures}" -gt 0 ]]; then
  echo "[hc-mem] completed with ${profile_failures} failed profile(s) (see report)" >&2
  # Smoke success + realistic pending is still a usable report; exit 0 so CI can publish docs.
  # Dedicated-host workflow should treat missing realistic numbers as a publish gate separately.
fi
