#!/usr/bin/env bash
# Reproducible HaplotypeCaller Peak-RSS profile: gatk-rs (Rust) vs pinned Java GATK 4.4.
#
# Output:
#   docs/perf/HC_MEMORY_PROFILE.md          — human-readable report (numbers + exact cmds)
#   docs/perf/hc_memory_profile_latest.json — machine-readable summary
#   docs/perf/runs/<stamp>/                — raw time logs + VCFs
#
# Requirements:
#   - cargo, samtools
#   - Java GATK via GATK_JAR, `gatk` on PATH, or Docker image from docs/GATK_PINNED.env
#   - macOS: /usr/bin/time -l ; Linux: /usr/bin/time -v (or gtime)
#
# Usage:
#   ./scripts/perf/run_hc_memory_profile.sh
#   JAVA_XMX=4g ./scripts/perf/run_hc_memory_profile.sh
#
# IMPORTANT: Do not invent README "X% less memory" claims without this report existing
# and being re-run after material engine changes.
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

# --- Fixture (tiny, checked-in; interval matches p4 assembly-region case) ---
REF="${HC_MEM_REF:-${repo_root}/parity/fixtures/reference.fa}"
BAM="${HC_MEM_BAM:-${repo_root}/parity/fixtures/sample.bam}"
INTERVAL="${HC_MEM_INTERVAL:-chr1:1-32}"
# Realistic pipeline heap for GATK HC shards (Broad WDLs commonly use 4g–6g).
JAVA_XMX="${JAVA_XMX:-4g}"
JAVA_XMS="${JAVA_XMS:-1g}"
JAVA_OPTS="${JAVA_OPTS:--Xms${JAVA_XMS} -Xmx${JAVA_XMX}}"

if [[ ! -f "${BAM}.bai" ]]; then
  samtools index "${BAM}"
fi

export CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-1}"
# Honour CARGO_TARGET_DIR when set (CI / sandbox); else default workspace target/.
target_dir="${CARGO_TARGET_DIR:-${repo_root}/target}"
echo "[hc-mem] building release gatk-cli (CARGO_TARGET_DIR=${target_dir})…"
(
  cd "${repo_root}"
  cargo build -p gatk-cli --release --locked
)
RUST_BIN="${target_dir}/release/gatk-rs"
if [[ ! -x "${RUST_BIN}" ]]; then
  echo "missing ${RUST_BIN}" >&2
  exit 1
fi

RUST_VCF="${run_dir}/rust.hc.vcf"
JAVA_VCF="${run_dir}/java.hc.vcf"
RUST_TIME="${run_dir}/rust.time.txt"
JAVA_TIME="${run_dir}/java.time.txt"
RUST_LOG="${run_dir}/rust.stdout.txt"
JAVA_LOG="${run_dir}/java.stdout.txt"

backend="$(giab_time_backend)"
echo "[hc-mem] time backend=${backend}"

# --- Rust native ---
echo "[hc-mem] Rust HaplotypeCaller…"
set +e
giab_run_timed "${backend}" "${RUST_TIME}" \
  "${RUST_BIN}" HaplotypeCaller \
    -R "${REF}" \
    -I "${BAM}" \
    -O "${RUST_VCF}" \
    -L "${INTERVAL}" \
  >"${RUST_LOG}" 2>&1
rust_exit=$?
set -e
# macOS time -l writes command stdout to .stdout; merge if present
if [[ -f "${RUST_TIME}.stdout" ]]; then
  cat "${RUST_TIME}.stdout" >>"${RUST_LOG}" || true
fi
if [[ "${rust_exit}" -ne 0 ]]; then
  echo "Rust HaplotypeCaller failed (exit=${rust_exit}). See ${RUST_LOG}" >&2
  tail -40 "${RUST_LOG}" >&2 || true
  exit "${rust_exit}"
fi

# --- Java GATK 4.4 (time measured *inside* container when using Docker) ---
echo "[hc-mem] Java GATK HaplotypeCaller (${GATK_DOCKER_IMAGE:-local})…"
java_cmd_file="${run_dir}/java.cmdline.txt"

run_java_timed_v2() {
  if [[ -n "${GATK_JAR:-}" && -f "${GATK_JAR}" ]]; then
    echo "java ${JAVA_OPTS} -jar ${GATK_JAR} HaplotypeCaller -R ${REF} -I ${BAM} -O ${JAVA_VCF} -L ${INTERVAL}" \
      | tee "${java_cmd_file}"
    giab_run_timed "${backend}" "${JAVA_TIME}" \
      java ${JAVA_OPTS} -jar "${GATK_JAR}" HaplotypeCaller \
        -R "${REF}" -I "${BAM}" -O "${JAVA_VCF}" -L "${INTERVAL}" \
      >"${JAVA_LOG}" 2>&1
    return $?
  fi
  if command -v gatk >/dev/null 2>&1; then
    echo "gatk --java-options '${JAVA_OPTS}' HaplotypeCaller -R ${REF} -I ${BAM} -O ${JAVA_VCF} -L ${INTERVAL}" \
      | tee "${java_cmd_file}"
    giab_run_timed "${backend}" "${JAVA_TIME}" \
      gatk --java-options "${JAVA_OPTS}" HaplotypeCaller \
        -R "${REF}" -I "${BAM}" -O "${JAVA_VCF}" -L "${INTERVAL}" \
      >"${JAVA_LOG}" 2>&1
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
    echo "  -R ${REF} -I ${BAM} -O ${JAVA_VCF} -L ${INTERVAL}"
  } | tee "${java_cmd_file}"

  # GATK 4.4 image has no GNU /usr/bin/time. Sample Peak RSS via /proc VmHWM while gatk runs.
  local inner="${run_dir}/java_inner.sh"
  cat >"${inner}" <<EOF
#!/bin/bash
set -euo pipefail
gatk --java-options '${JAVA_OPTS}' HaplotypeCaller \
  -R '${REF}' -I '${BAM}' -O '${JAVA_VCF}' -L '${INTERVAL}' \
  >'${JAVA_LOG}' 2>&1 &
gp=\$!
start_ns=\$(date +%s%N)
peak_kb=0
while kill -0 "\$gp" 2>/dev/null; do
  for status in /proc/[0-9]*/status; do
    [[ -r "\$status" ]] || continue
    # Only count java / gatk wrapper (not the sampler shell).
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
# GNU-time-compatible fields for giab_parse_time_log (wall as m:ss).
{
  echo "Elapsed (wall clock) time (h:mm:ss or m:ss): 0:\${wall}"
  echo "Maximum resident set size (kbytes): \${peak_kb}"
  echo "Command being timed: gatk HaplotypeCaller (Linux /proc VmHWM sampler; image lacks GNU time)"
} >'${JAVA_TIME}'
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
run_java_timed_v2
java_exit=$?
set -e
if [[ "${java_exit}" -ne 0 ]]; then
  echo "Java HaplotypeCaller failed (exit=${java_exit}). See ${JAVA_LOG}" >&2
  tail -60 "${JAVA_LOG}" >&2 || true
  exit "${java_exit}"
fi

rust_json="$(giab_parse_time_log "${RUST_TIME}")"
java_json="$(giab_parse_time_log "${JAVA_TIME}")"

rustc_version="$(rustc --version 2>/dev/null || echo unknown)"
cargo_version="$(cargo --version 2>/dev/null || echo unknown)"
uname_s="$(uname -srm)"
rust_sha="$(git -C "${repo_root}" rev-parse --short HEAD 2>/dev/null || echo unknown)"

python3 - "${run_dir}" "${out_root}" "${stamp}" \
  "${REF}" "${BAM}" "${INTERVAL}" \
  "${JAVA_OPTS}" "${GATK_DOCKER_IMAGE:-}" "${GATK_PINNED_SHA:-}" \
  "${rust_json}" "${java_json}" \
  "${rustc_version}" "${cargo_version}" "${uname_s}" "${rust_sha}" \
  "${RUST_BIN}" <<'PY'
import json, os, sys, pathlib, datetime

(
    run_dir, out_root, stamp,
    ref, bam, interval,
    java_opts, docker_image, gatk_sha,
    rust_json, java_json,
    rustc_version, cargo_version, uname_s, rust_sha,
    rust_bin,
) = sys.argv[1:]

rust = json.loads(rust_json)
java = json.loads(java_json)

def fmt_mib(kb):
    if kb is None:
        return "n/a"
    return f"{kb / 1024.0:.2f} MiB ({kb} KiB)"

summary = {
    "stamp_utc": stamp,
    "host": uname_s,
    "fixture": {
        "reference": ref,
        "bam": bam,
        "interval": interval,
        "note": "Checked-in p4 smoke fixture (chr1 32bp). Peak RSS on this tiny region is dominated by process/runtime overhead; do not market as genome-wide memory savings.",
    },
    "rust": {
        "binary": rust_bin,
        "git_sha": rust_sha,
        "build": "cargo build -p gatk-cli --release --locked",
        "rustc": rustc_version,
        "cargo": cargo_version,
        "wall_sec": rust.get("wall_sec"),
        "max_rss_kb": rust.get("max_rss_kb"),
        "max_rss_mib": None if rust.get("max_rss_kb") is None else rust["max_rss_kb"] / 1024.0,
        "cmdline": f"{rust_bin} HaplotypeCaller -R {ref} -I {bam} -O <run>/rust.hc.vcf -L {interval}",
    },
    "java": {
        "gatk_version_pin": "4.4.0.0",
        "gatk_sha": gatk_sha,
        "docker_image": docker_image or None,
        "java_options": java_opts,
        "wall_sec": java.get("wall_sec"),
        "max_rss_kb": java.get("max_rss_kb"),
        "max_rss_mib": None if java.get("max_rss_kb") is None else java["max_rss_kb"] / 1024.0,
        "note": "Peak RSS = max /proc VmHWM (kB) of java|gatk inside Docker when the image lacks GNU time; else GNU time -v on the JVM process.",
    },
}

# Optional ratio — only if both present; never invent marketing %.
rk, jk = rust.get("max_rss_kb"), java.get("max_rss_kb")
if rk and jk and jk > 0:
    summary["ratio"] = {
        "java_over_rust_peak_rss": jk / rk,
        "rust_fraction_of_java_peak_rss": rk / jk,
        "delta_kb": jk - rk,
    }

pathlib.Path(out_root).mkdir(parents=True, exist_ok=True)
json_path = pathlib.Path(out_root) / "hc_memory_profile_latest.json"
json_path.write_text(json.dumps(summary, indent=2) + "\n")
(pathlib.Path(run_dir) / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")

ratio_block = ""
if "ratio" in summary:
    r = summary["ratio"]
    ratio_block = f"""
| Java / Rust Peak-RSS | {r['java_over_rust_peak_rss']:.2f}× |
| Rust as fraction of Java Peak-RSS | {100.0 * r['rust_fraction_of_java_peak_rss']:.1f}% |
| Absolute delta (Java − Rust) | {r['delta_kb']/1024.0:.2f} MiB |
"""

md = f"""# HaplotypeCaller memory profile (reproducible)

**Generated (UTC):** `{stamp}`  
**Host:** `{uname_s}`  
**Runner script:** [`scripts/perf/run_hc_memory_profile.sh`](../../scripts/perf/run_hc_memory_profile.sh)  
**Raw run directory:** `docs/perf/runs/{stamp}/`

> **Scope warning:** This profile uses the checked-in p4 smoke fixture
> (`parity/fixtures/sample.bam` + `reference.fa`, interval `{interval}`).
> Absolute Peak-RSS is dominated by runtime/JVM fixed costs on such a tiny
> window. **Do not** advertise these numbers as genome-wide “X% less memory”
> without re-measuring on a realistic GIAB shard.

## Peak-RSS (side by side)

| Engine | Peak RSS | Wall time |
|--------|----------|-----------|
| **gatk-rs** (Rust release) | **{fmt_mib(rk)}** | {rust.get('wall_sec')} s |
| **Java GATK 4.4.0.0** | **{fmt_mib(jk)}** | {java.get('wall_sec')} s |
{ratio_block}

## Exact commands

### Rust

```bash
cargo build -p gatk-cli --release --locked
# rustc: {rustc_version}
# cargo: {cargo_version}
# git: {rust_sha}
{rust_bin} HaplotypeCaller \\
  -R {ref} \\
  -I {bam} \\
  -O docs/perf/runs/{stamp}/rust.hc.vcf \\
  -L {interval}
```

Time capture (this host): see `docs/perf/runs/{stamp}/rust.time.txt`
(macOS `/usr/bin/time -l` or GNU `/usr/bin/time -v`).

### Java GATK 4.4

- Pin: `GATK_PINNED_SHA={gatk_sha}` (`docs/GATK_PINNED.env`)
- Image: `{docker_image or 'local gatk / GATK_JAR'}`
- JVM options (pipeline-realistic): `{java_opts}`

```bash
# Re-run via the harness (preferred):
./scripts/perf/run_hc_memory_profile.sh
# Or see docs/perf/runs/{stamp}/java.cmdline.txt for the exact docker/java line used.
```

Time capture: `docs/perf/runs/{stamp}/java.time.txt`  
When Docker is used, Peak-RSS is sampled from `/proc/*/status` **VmHWM**
for `java`/`gatk` **inside** the Linux container (the Broad 4.4 image has no
GNU `/usr/bin/time`). Host `time docker …` is never used for RSS.

## Optional deeper profiling

- **macOS Instruments:**  
  `xcrun xctrace record --template 'Allocations' --output docs/perf/runs/{stamp}/rust.allocations.trace --launch -- {rust_bin} HaplotypeCaller -R {ref} -I {bam} -O /tmp/rust.hc.vcf -L {interval}`
- **Linux heaptrack (Docker):**  
  mount the release binary and fixture into a heaptrack image; keep Peak-RSS
  from this script as the primary comparable number.

## Re-run

```bash
./scripts/perf/run_hc_memory_profile.sh
# optional overrides:
#   JAVA_XMX=4g JAVA_XMS=1g HC_MEM_INTERVAL=chr1:1-32 ./scripts/perf/run_hc_memory_profile.sh
```
"""
md_path = pathlib.Path(out_root) / "HC_MEMORY_PROFILE.md"
md_path.write_text(md)
print(f"Wrote {md_path}")
print(f"Wrote {json_path}")
print(f"Rust Peak-RSS KiB={rk}  Java Peak-RSS KiB={jk}")
PY

echo "[hc-mem] done → docs/perf/HC_MEMORY_PROFILE.md"
