#!/usr/bin/env bash
# L8 FORMAT gate — hard GT + soft AD/DP among GT-matched dense truth TPs.
#
# Soft contract (among GT-matched): AD L1≤2, |DP|≤2.
# Exact PL remains waived (W-L7-FORMAT); tracked separately if L8_FORMAT_INCLUDE_PL=1.
#
# Usage:
#   export P12_REFERENCE="$PWD/parity/realworld/assets/hs37d5.simple.fa"
#   ./scripts/parity/run_l8_format_gate.sh
#
# Optional:
#   L8_FORMAT_INTERVAL=20:10000000-10050000
#   L8_FORMAT_JAVA_VCF / L8_FORMAT_RUST_VCF
#   L8_MAX_HARD_MISMATCH_RATE=0.15
#   L8_MAX_SOFT_MISMATCH_RATE=0.30   # AD/DP among GT-matched
#   L8_FORMAT_INCLUDE_PL=1          # also gate PL with --pl-tol
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"

assets="${repo_root}/parity/realworld/assets"
report_dir="${L8_FORMAT_REPORT_DIR:-${repo_root}/parity/reports/hc-full-parity-j6-dense}"
interval="${L8_FORMAT_INTERVAL:-20:10000000-10050000}"
java_vcf="${L8_FORMAT_JAVA_VCF:-${report_dir}/p12_dense_giab_window.java.vcf}"
rust_vcf="${L8_FORMAT_RUST_VCF:-${report_dir}/p12_dense_giab_window.rust.vcf}"
truth_vcf="${L8_FORMAT_TRUTH_VCF:-${assets}/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz}"
regions_bed="${L8_FORMAT_REGIONS_BED:-${assets}/HG001_GRCh37_1_22_v4.2.1_benchmark.bed}"

hard_max="${L8_MAX_HARD_MISMATCH_RATE:-0.15}"
soft_max="${L8_MAX_SOFT_MISMATCH_RATE:-0.30}"
soft_keys="AD,DP"
if [[ "${L8_FORMAT_INCLUDE_PL:-0}" == "1" ]]; then
  soft_keys="AD,DP,PL"
fi

stamp="$(date -u +%Y%m%dT%H%M%SZ)"
json_out="${report_dir}/l8_format_gate_${stamp}.json"
md_out="${report_dir}/l8_format_gate_${stamp}.md"
canonical_json="${report_dir}/l8_format_gate.json"
canonical_md="${report_dir}/l8_format_gate.md"

if [[ ! -f "${java_vcf}" || ! -f "${rust_vcf}" ]]; then
  echo "[l8-format] missing VCFs; run J6_DENSE=1 ./scripts/parity/run_hc_full_parity_j6_truth.sh first" >&2
  echo "  java=${java_vcf}" >&2
  echo "  rust=${rust_vcf}" >&2
  exit 1
fi

echo "=== L8 FORMAT gate ${stamp} ==="
echo "interval=${interval}"
echo "soft_keys=${soft_keys} soft_max=${soft_max} hard_max=${hard_max}"

python3 "${repo_root}/scripts/parity/l7_dense_format_spotcheck.py" \
  --java-vcf "${java_vcf}" \
  --rust-vcf "${rust_vcf}" \
  --truth-vcf "${truth_vcf}" \
  --regions-bed "${regions_bed}" \
  --eval-interval "${interval}" \
  --max-sites "${L8_FORMAT_MAX_SITES:-119}" \
  --hard-keys GT \
  --soft-keys "${soft_keys}" \
  --ad-l1-tol "${L8_AD_L1_TOL:-2}" \
  --dp-tol "${L8_DP_TOL:-2}" \
  --pl-tol "${L8_PL_TOL:-5}" \
  --max-hard-mismatch-rate "${hard_max}" \
  --max-soft-mismatch-rate "${soft_max}" \
  --json-out "${json_out}" \
  --md-out "${md_out}"

cp -f "${json_out}" "${canonical_json}"
cp -f "${md_out}" "${canonical_md}"

echo "=== L8 FORMAT gate PASS ==="
echo "json: ${canonical_json}"
echo "md: ${canonical_md}"
