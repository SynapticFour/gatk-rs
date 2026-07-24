#!/usr/bin/env bash
# Fair HC comparison on the dedicated quiet host.
#
# Configs (each × PERF_REPEATS, default 5):
#   rust_logless_scalar  — gatk-rs --pair-hmm LOGLESS_HMM
#   rust_simd            — gatk-rs --pair-hmm AVX
#   java_fastest_available — Java --pair-hmm-implementation FASTEST_AVAILABLE
#   java_logless_caching — Java --pair-hmm-implementation LOGLESS_CACHING
#
# Regions (nested on GIAB dense window when staged):
#   small / medium / large
#
# Metrics: wall, user, sys, Peak-RSS; optional RAPL energy via perf stat.
#
# Env:
#   PERF_REPEATS=5
#   PERF_THREADS=1
#   PERF_CPU_LIST=0-3
#   PERF_SKIP_ENERGY=0
#   PERF_STAGE_GIAB=1   — stage NA12878 GIAB window BAM if missing
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"
# shellcheck source=../parity/lib_pinned_gatk.sh
source "${repo_root}/scripts/parity/lib_pinned_gatk.sh"
# shellcheck source=../parity/giab/lib_giab.sh
source "${repo_root}/scripts/parity/giab/lib_giab.sh"

REPEATS="${PERF_REPEATS:-5}"
THREADS="${PERF_THREADS:-1}"
JAVA_XMX="${JAVA_XMX:-4g}"
JAVA_XMS="${JAVA_XMS:-1g}"
JAVA_OPTS="${JAVA_OPTS:--Xms${JAVA_XMS} -Xmx${JAVA_XMX}}"
stamp="$(date -u +%Y%m%dT%H%M%SZ)"
out_root="${repo_root}/docs/perf"
run_dir="${PERF_RUN_DIR:-${out_root}/runs/fair_${stamp}}"
raw_dir="${run_dir}/raw"
mkdir -p "${raw_dir}" "${out_root}"

# Pin the *measured* process (not a bash function — taskset cannot exec shell funcs).
timed_run() {
  local backend="$1" time_log="$2"
  shift 2
  if [[ -n "${PERF_CPU_LIST:-}" ]] && command -v taskset >/dev/null 2>&1; then
    giab_run_timed "${backend}" "${time_log}" taskset -c "${PERF_CPU_LIST}" "$@"
  else
    giab_run_timed "${backend}" "${time_log}" "$@"
  fi
}

echo "[fair-hc] capturing host specs…"
export PERF_RUN_DIR="${run_dir}"
export PERF_RUNNER_LABEL="${PERF_RUNNER_LABEL:-gatk-rs-benchmark}"
"${script_dir}/capture_host_specs.sh"

# --- Stage inputs (prefer dense GIAB window for meaningful region sizes) ---
REF="${PERF_REF:-${repo_root}/parity/realworld/assets/hs37d5.simple.fa}"
BAM_DIR="${repo_root}/parity/realworld/na12878_giab_window_b37"
BAM="${PERF_BAM:-${BAM_DIR}/NA12878_giab_window.b37.bam}"

if [[ ! -f "${REF}" ]]; then
  echo "[fair-hc] staging reference…"
  "${repo_root}/scripts/parity/realworld/03_stage_reference_and_truth.sh"
fi
if [[ ! -f "${BAM}" ]]; then
  if [[ "${PERF_STAGE_GIAB:-1}" == "1" ]]; then
    echo "[fair-hc] staging GIAB dense window BAM…"
    "${repo_root}/scripts/parity/realworld/04_stage_na12878_giab_window_bam.sh"
  else
    echo "[fair-hc] ERROR: missing BAM ${BAM}" >&2
    exit 1
  fi
fi
if [[ ! -f "${BAM}.bai" && ! -f "${BAM%.*}.bai" ]]; then
  samtools index "${BAM}"
fi

# Nested intervals on the default GIAB window (20:10000000-10050000).
INTERVAL_SMALL="${PERF_INTERVAL_SMALL:-20:10000000-10005000}"
INTERVAL_MEDIUM="${PERF_INTERVAL_MEDIUM:-20:10000000-10020000}"
INTERVAL_LARGE="${PERF_INTERVAL_LARGE:-20:10000000-10050000}"

declare -a REGION_IDS=(small medium large)
declare -a REGION_INTERVALS=("${INTERVAL_SMALL}" "${INTERVAL_MEDIUM}" "${INTERVAL_LARGE}")

# Build Rust binary once.
export CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-4}"
target_dir="${CARGO_TARGET_DIR:-${repo_root}/target}"
echo "[fair-hc] building release gatk-cli…"
(
  cd "${repo_root}"
  cargo build -p gatk-cli --release --locked
)
RUST_BIN="${target_dir}/release/gatk-rs"
[[ -x "${RUST_BIN}" ]] || { echo "missing ${RUST_BIN}" >&2; exit 1; }

# Verify Java native path once before timing (fail closed for published runs).
if [[ "${PERF_SKIP_JAVA_AVX_VERIFY:-0}" != "1" ]]; then
  echo "[fair-hc] verifying Java FASTEST_AVAILABLE uses native PairHMM…"
  HC_MEM_REF="${REF}" HC_MEM_BAM="${BAM}" HC_MEM_INTERVAL="${INTERVAL_SMALL}" \
    "${script_dir}/verify_java_avx_pairhmm.sh"
fi

backend="$(giab_time_backend)"
echo "[fair-hc] time_backend=${backend} repeats=${REPEATS} threads=${THREADS}"

# Optional energy wrapper (RAPL). Soft-fail if unavailable / needs privileges.
energy_prefix() {
  # prints nothing; sets global _energy_cmd array via nameref-ish pattern
  _use_perf=0
  if [[ "${PERF_SKIP_ENERGY:-0}" == "1" ]]; then
    return 0
  fi
  if ! command -v perf >/dev/null 2>&1; then
    return 0
  fi
  # Probe: does energy-pkg exist?
  if perf list 2>/dev/null | grep -q 'power/energy-pkg'; then
    _use_perf=1
  fi
}

run_one() {
  local config_id="$1" label="$2" engine="$3" pair_hmm="$4" region_size="$5" interval="$6" rep="$7"
  local cell_dir="${raw_dir}/${region_size}/${config_id}/rep$(printf '%02d' "${rep}")"
  mkdir -p "${cell_dir}"
  local out_vcf="${cell_dir}/out.vcf"
  local time_log="${cell_dir}/time.txt"
  local stdout_log="${cell_dir}/stdout.txt"
  local perf_log="${cell_dir}/perf.txt"
  local metrics_json="${cell_dir}/metrics.json"
  local exit_code=0

  echo "[fair-hc] ${region_size} ${config_id} rep=${rep}/${REPEATS}"

  # Drop caches optionally between runs? Too invasive for cloud VMs — skip.

  local -a cmd=()
  if [[ "${engine}" == "rust" ]]; then
    cmd=(
      "${RUST_BIN}" HaplotypeCaller
      -R "${REF}" -I "${BAM}" -O "${out_vcf}" -L "${interval}"
      --threads "${THREADS}"
      --pair-hmm "${pair_hmm}"
    )
  else
    # Java via Docker (pinned image) — time wraps docker so overhead is equal across Java configs.
    cmd=(docker run --rm)
    if [[ "$(uname -s)" == "Darwin" ]]; then
      cmd+=(--platform "${GATK_DOCKER_PLATFORM:-linux/amd64}")
    fi
    cmd+=(
      -v "${repo_root}:${repo_root}"
      -w "${repo_root}"
      "${GATK_DOCKER_IMAGE}"
      gatk --java-options "${JAVA_OPTS}" HaplotypeCaller
      -R "${REF}" -I "${BAM}" -O "${out_vcf}" -L "${interval}"
      --pair-hmm-implementation "${pair_hmm}"
      --native-pair-hmm-threads "${THREADS}"
    )
  fi

  set +e
  if [[ "${_use_perf:-0}" == "1" ]]; then
    timed_run "${backend}" "${time_log}" \
      perf stat -o "${perf_log}" -e power/energy-pkg/,power/energy-cores/ -- \
      "${cmd[@]}" >"${stdout_log}" 2>&1
    exit_code=$?
  else
    timed_run "${backend}" "${time_log}" \
      "${cmd[@]}" >"${stdout_log}" 2>&1
    exit_code=$?
    : >"${perf_log}"
  fi
  set -e

  if [[ -f "${time_log}.stdout" ]]; then
    cat "${time_log}.stdout" >>"${stdout_log}" || true
  fi

  local parsed
  parsed="$(python3 "${script_dir}/parse_time_metrics.py" "${time_log}" --perf-log "${perf_log}")"

  python3 - "${metrics_json}" "${parsed}" <<PY
import json, sys
path, parsed_s = sys.argv[1], sys.argv[2]
m = json.loads(parsed_s)
m.update({
    "config_id": "${config_id}",
    "label": "${label}",
    "engine": "${engine}",
    "pair_hmm": "${pair_hmm}",
    "region_size": "${region_size}",
    "interval": "${interval}",
    "repeat": int("${rep}"),
    "ok": ${exit_code} == 0,
    "exit_code": int("${exit_code}"),
    "threads": int("${THREADS}"),
    "gatk_docker": "${GATK_DOCKER_IMAGE}",
    "gatk_pinned_sha": "${GATK_PINNED_SHA}",
})
open(path, "w", encoding="utf-8").write(json.dumps(m, indent=2) + "\n")
PY

  if [[ "${exit_code}" -ne 0 ]]; then
    echo "[fair-hc] WARN: ${config_id}/${region_size} rep=${rep} failed (exit=${exit_code})" >&2
    tail -30 "${stdout_log}" >&2 || true
  fi
}

energy_prefix

# Config table
# id|label|engine|pair_hmm
CONFIGS=(
  "rust_logless_scalar|gatk-rs LOGLESS_HMM scalar|rust|LOGLESS_HMM"
  "rust_simd|gatk-rs AVX/SIMD|rust|AVX"
  "java_fastest_available|Java GATK FASTEST_AVAILABLE|java|FASTEST_AVAILABLE"
  "java_logless_caching|Java GATK LOGLESS_CACHING|java|LOGLESS_CACHING"
)

# Warm-up discarded (one short run) to reduce first-JIT / page-cache shock.
if [[ "${PERF_SKIP_WARMUP:-0}" != "1" ]]; then
  echo "[fair-hc] warmup (discarded)…"
  if [[ -n "${PERF_CPU_LIST:-}" ]] && command -v taskset >/dev/null 2>&1; then
    taskset -c "${PERF_CPU_LIST}" "${RUST_BIN}" HaplotypeCaller \
      -R "${REF}" -I "${BAM}" -O "${run_dir}/warmup.vcf" -L "${INTERVAL_SMALL}" \
      --threads "${THREADS}" --pair-hmm LOGLESS_HMM \
      >"${run_dir}/warmup.log" 2>&1 || true
  else
    "${RUST_BIN}" HaplotypeCaller \
      -R "${REF}" -I "${BAM}" -O "${run_dir}/warmup.vcf" -L "${INTERVAL_SMALL}" \
      --threads "${THREADS}" --pair-hmm LOGLESS_HMM \
      >"${run_dir}/warmup.log" 2>&1 || true
  fi
fi

for idx in "${!REGION_IDS[@]}"; do
  region="${REGION_IDS[$idx]}"
  interval="${REGION_INTERVALS[$idx]}"
  for cfg_line in "${CONFIGS[@]}"; do
    IFS='|' read -r cid clabel ceng cph <<<"${cfg_line}"
    for ((rep = 1; rep <= REPEATS; rep++)); do
      run_one "${cid}" "${clabel}" "${ceng}" "${cph}" "${region}" "${interval}" "${rep}"
    done
  done
done

meta_json="${run_dir}/meta.json"
python3 - "${meta_json}" <<PY
import json, os
open("${meta_json}", "w").write(json.dumps({
    "repeats": int("${REPEATS}"),
    "threads": int("${THREADS}"),
    "commit_sha": os.environ.get("GITHUB_SHA", "local"),
    "workflow_run_url": os.environ.get("PERF_WORKFLOW_RUN_URL", ""),
    "ref": "${REF}",
    "bam": "${BAM}",
    "intervals": {
        "small": "${INTERVAL_SMALL}",
        "medium": "${INTERVAL_MEDIUM}",
        "large": "${INTERVAL_LARGE}",
    },
    "java_xmx": "${JAVA_XMX}",
    "gatk_docker": "${GATK_DOCKER_IMAGE}",
    "gatk_pinned_sha": "${GATK_PINNED_SHA}",
    "primary_java_baseline": "FASTEST_AVAILABLE",
}, indent=2) + "\n")
PY

echo "[fair-hc] aggregating…"
python3 "${script_dir}/aggregate_fair_comparison.py" \
  --raw-dir "${raw_dir}" \
  --host-json "${run_dir}/host_specs.json" \
  --meta-json "${meta_json}" \
  --out-json "${out_root}/fair_hc_comparison_latest.json" \
  --out-md "${out_root}/FAIR_HC_COMPARISON.md"

cp -f "${out_root}/fair_hc_comparison_latest.json" "${run_dir}/fair_hc_comparison.json"
cp -f "${out_root}/FAIR_HC_COMPARISON.md" "${run_dir}/FAIR_HC_COMPARISON.md"

# Pointer
cat >"${out_root}/DEDICATED_BENCH_LATEST.md" <<EOF
# Latest dedicated-host benchmark pointer

**Fair HC comparison:** [\`FAIR_HC_COMPARISON.md\`](FAIR_HC_COMPARISON.md)  
**JSON:** [\`fair_hc_comparison_latest.json\`](fair_hc_comparison_latest.json)  
**Run dir:** [\`runs/fair_${stamp}/\`](runs/fair_${stamp}/)  
**Host:** [\`HOST_SPECS.md\`](HOST_SPECS.md)

Primary Java baseline for speedups: **\`FASTEST_AVAILABLE\`** (not \`LOGLESS_CACHING\`).
EOF

echo "[fair-hc] done → ${run_dir}"
echo "[fair-hc] report → ${out_root}/FAIR_HC_COMPARISON.md"
