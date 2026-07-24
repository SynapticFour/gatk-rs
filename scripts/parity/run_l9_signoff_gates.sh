#!/usr/bin/env bash
# L9 production sign-off gate battery (FORMAT + F1 on dense slices; optional P12).
#
# Assumes release HC binary already built and dense report VCFs exist (or regenerate
# with P12_SKIP_JAVA=1 when java VCFs are cached).
#
# Usage (M4-safe):
#   source scripts/parity/m4_disk_guard.sh && m4_require_free_gb 12
#   export CARGO_TARGET_DIR="$PWD/target" CARGO_BUILD_JOBS=1 RAYON_NUM_THREADS=2
#   export P12_REFERENCE="$PWD/parity/realworld/assets/hs37d5.simple.fa"
#   ./scripts/parity/run_l9_signoff_gates.sh
#
# Env:
#   L9_RUN_P12=1           also run ./scripts/parity/run_p12_l3_signoff.sh
#   L9_REGEN_RUST=1        re-run dense HC (rust) before FORMAT/F1 (java cached)
#   L9_INCLUDE_PL_INFO=1   informational soft PL stats (does not fail L9)
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"

if [[ -f "${repo_root}/scripts/parity/m4_disk_guard.sh" ]]; then
  # shellcheck disable=SC1091
  source "${repo_root}/scripts/parity/m4_disk_guard.sh"
  m4_require_free_gb "${L9_MIN_FREE_GB:-10}"
fi

export P12_REFERENCE="${P12_REFERENCE:-${repo_root}/parity/realworld/assets/hs37d5.simple.fa}"
export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-${repo_root}/target}"
export CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-1}"
export RAYON_NUM_THREADS="${RAYON_NUM_THREADS:-2}"

stamp="$(date -u +%Y%m%dT%H%M%SZ)"
log="${repo_root}/parity/reports/l9_signoff_gates_${stamp}.log"
mkdir -p "${repo_root}/parity/reports"
exec > >(tee -a "${log}") 2>&1

echo "=== L9 sign-off gates ${stamp} ==="
echo "P12_REFERENCE=${P12_REFERENCE}"
echo "CARGO_TARGET_DIR=${CARGO_TARGET_DIR}"

run_slice_format() {
  local name="$1"
  local interval="$2"
  local report_dir="$3"
  local max_sites="${4:-200}"
  echo ""
  echo "=== FORMAT ${name} (${interval}) ==="
  L8_FORMAT_INTERVAL="${interval}" \
  L8_FORMAT_REPORT_DIR="${report_dir}" \
  L8_FORMAT_MAX_SITES="${max_sites}" \
  L8_FORMAT_INCLUDE_PL="${L9_INCLUDE_PL_INFO:-0}" \
    ./scripts/parity/run_l8_format_gate.sh
}

run_slice_f1() {
  local name="$1"
  local report_dir="$2"
  local thresholds="$3"
  local interval="$4"
  echo ""
  echo "=== F1 ${name} ==="
  python3 "${repo_root}/scripts/parity/p13_truth_eval.py" \
    --java-vcf "${report_dir}/p12_dense_giab_window.java.vcf" \
    --rust-vcf "${report_dir}/p12_dense_giab_window.rust.vcf" \
    --truth-vcf "${repo_root}/parity/realworld/assets/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz" \
    --regions-bed "${repo_root}/parity/realworld/assets/HG001_GRCh37_1_22_v4.2.1_benchmark.bed" \
    --eval-interval "${interval}" \
    --thresholds-json "${thresholds}" \
    --json-out "${report_dir}/l9_f1_gate.json" \
    --md-out "${report_dir}/l9_f1_gate.md"
}

if [[ "${L9_REGEN_RUST:-0}" == "1" ]]; then
  echo ""
  echo "=== regen dense rust VCFs (java cached) ==="
  export P12_SKIP_JAVA=1
  J6_DENSE=1 ./scripts/parity/run_hc_full_parity_j6_truth.sh
  ./scripts/parity/run_hc_full_parity_j6_dense_chr21.sh
  ./scripts/parity/run_hc_full_parity_j6_dense_holdout.sh
fi

# chr20 primary dense
run_slice_format "chr20" "20:10000000-10050000" \
  "${repo_root}/parity/reports/hc-full-parity-j6-dense" 119
run_slice_f1 "chr20" \
  "${repo_root}/parity/reports/hc-full-parity-j6-dense" \
  "${repo_root}/parity/fixtures/hc-full-parity/j6/thresholds_dense.json" \
  "20:10000000-10050000"

# chr21
run_slice_format "chr21" "21:41200001-41250000" \
  "${repo_root}/parity/reports/hc-full-parity-j6-dense-chr21" 200
run_slice_f1 "chr21" \
  "${repo_root}/parity/reports/hc-full-parity-j6-dense-chr21" \
  "${repo_root}/parity/fixtures/hc-full-parity/j6/thresholds_dense_chr21.json" \
  "21:41200001-41250000"

# holdout (W-L9-HOLDOUT-INDEL on stratified indel floors)
run_slice_format "holdout" "20:15000000-15050000" \
  "${repo_root}/parity/reports/hc-full-parity-j6-dense-holdout" 200
run_slice_f1 "holdout" \
  "${repo_root}/parity/reports/hc-full-parity-j6-dense-holdout" \
  "${repo_root}/parity/fixtures/hc-full-parity/j6/thresholds_dense_holdout.json" \
  "20:15000000-15050000"

if [[ "${L9_RUN_P12:-0}" == "1" ]]; then
  echo ""
  echo "=== P12 L3 ==="
  ./scripts/parity/run_p12_l3_signoff.sh
fi

echo ""
echo "=== L9 sign-off gates PASS ==="
echo "log: ${log}"
echo "See docs/CLAIM_MATRIX.md (W-L7-FORMAT, W-L9-HOLDOUT-INDEL)."
