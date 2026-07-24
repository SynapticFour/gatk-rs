#!/usr/bin/env bash
# L6 gate: NA12878 scale (P12 interval) + GIAB truth stratification (P13).
#
# Usage:
#   export P12_REFERENCE="$PWD/parity/realworld/assets/hs37d5.simple.fa"
#   ./scripts/parity/run_hc_full_parity_j6_truth.sh
#
# Optional:
#   J6_INTERVAL=2:92000000-92400000   # default GIAB-overlapping window (includes P12 spine)
#   J6_SKIP_HC=1                      # reuse existing VCFs in parity/reports/
#   J6_STRICT=0                       # report only, do not fail on gate
#   J6_DENSE=1                        # R3: stage dense GIAB window BAM + non-vacuous probe
#   J6_REPORT_DIR=...                 # L7: override report tree (second dense slices)
#   J6_DENSE_OUT_DIR=...              # L7: override dense BAM staging directory
#   J6_THRESHOLDS=...                 # override thresholds JSON
#
# Reports: parity/reports/hc-full-parity-j6/  (or hc-full-parity-j6-dense/ when J6_DENSE=1)
# Second slice helper: ./scripts/parity/run_hc_full_parity_j6_dense_chr21.sh
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"
cd "${repo_root}"

assets_dir="${J6_ASSETS_DIR:-${repo_root}/parity/realworld/assets}"
dense="${J6_DENSE:-0}"
if [[ "${dense}" == "1" ]]; then
  thresholds="${J6_THRESHOLDS:-${repo_root}/parity/fixtures/hc-full-parity/j6/thresholds_dense.json}"
  # L7-A3: override with J6_REPORT_DIR for additional non-chr20 dense slices.
  report_dir="${J6_REPORT_DIR:-${repo_root}/parity/reports/hc-full-parity-j6-dense}"
  eval_interval="${J6_INTERVAL:-${P12_INTERVAL:-20:10000000-10050000}}"
else
  thresholds="${J6_THRESHOLDS:-${repo_root}/parity/fixtures/hc-full-parity/j6/thresholds.json}"
  report_dir="${J6_REPORT_DIR:-${repo_root}/parity/reports/hc-full-parity-j6}"
  eval_interval="${J6_INTERVAL:-${P12_INTERVAL:-2:92000000-92400000}}"
fi
mkdir -p "${report_dir}"

strict="${J6_STRICT:-1}"

export P12_REFERENCE="${P12_REFERENCE:-${assets_dir}/hs37d5.simple.fa}"
export P13_TRUTH_VCF="${P13_TRUTH_VCF:-${assets_dir}/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz}"
export P13_REGIONS_BED="${P13_REGIONS_BED:-${assets_dir}/HG001_GRCh37_1_22_v4.2.1_benchmark.bed}"
export P12_INTERVAL="${eval_interval}"
export P13_EVAL_INTERVAL="${P13_EVAL_INTERVAL:-${eval_interval}}"
export PARITY_CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"
export P12_CARGO_RELEASE="${P12_CARGO_RELEASE:-1}"

stamp="$(date -u +%Y%m%dT%H%M%SZ)"
log="${report_dir}/j6_truth_${stamp}.log"
exec > >(tee -a "${log}") 2>&1

echo "=== L6 hc-full-parity j6 truth ${stamp} ==="
echo "P12_REFERENCE=${P12_REFERENCE}"
echo "P13_TRUTH_VCF=${P13_TRUTH_VCF}"
echo "eval_interval=${eval_interval}"
echo "thresholds=${thresholds}"

if [[ ! -f "${P12_REFERENCE}" ]]; then
  echo "[j6-truth] reference missing: ${P12_REFERENCE}" >&2
  echo "[j6-truth] stage with: ./scripts/parity/realworld/03_stage_reference_and_truth.sh" >&2
  exit 1
fi

if [[ ! -f "${P13_TRUTH_VCF}" ]]; then
  echo "[j6-truth] truth VCF missing: ${P13_TRUTH_VCF}" >&2
  echo "[j6-truth] stage with: ./scripts/parity/realworld/03_stage_reference_and_truth.sh" >&2
  exit 1
fi

if [[ "${dense}" == "1" ]]; then
  p12_json="${report_dir}/p12_dense_giab_window.json"
  p13_json="${report_dir}/p13_dense_truth_eval.json"
  export P12_JAVA_VCF="${report_dir}/p12_dense_giab_window.java.vcf"
  export P12_RUST_VCF="${report_dir}/p12_dense_giab_window.rust.vcf"
  export P12_JSON_OUT="${p12_json}"
  export P12_MD_OUT="${report_dir}/p12_dense_giab_window.md"
  export P13_JAVA_VCF="${P12_JAVA_VCF}"
  export P13_RUST_VCF="${P12_RUST_VCF}"
else
  p12_json="${repo_root}/parity/reports/p12_realworld_na12878_20k.json"
  p13_json="${repo_root}/parity/reports/p13_realworld_truth_eval.json"
fi
# Staging may rewrite p12 JSON to a "reference missing" stub — preserve when reusing HC.
p12_json_backup=""
if [[ "${J6_SKIP_HC:-0}" == "1" && -f "${p12_json}" ]]; then
  p12_json_backup="$(mktemp)"
  cp -f "${p12_json}" "${p12_json_backup}"
fi

echo ""
if [[ "${dense}" == "1" ]]; then
  echo "=== stage dense GIAB window BAM (R3) ==="
  ./scripts/parity/realworld/04_stage_na12878_giab_window_bam.sh
  dense_bam="${J6_DENSE_OUT_DIR:-${repo_root}/parity/realworld/na12878_giab_window_b37}/NA12878_giab_window.b37.bam"
  # Dense mode must not inherit a sparse P12_BAM from a prior L6 run in the same shell.
  # Override via J6_DENSE_BAM only; plain P12_BAM is for non-dense L6.
  export P12_BAM="${J6_DENSE_BAM:-${dense_bam}}"
  export P12_BAI="${P12_BAI:-${P12_BAM%.bam}.bai}"
  echo ""
  echo "=== vacuity probe (require non-vacuous GIAB depth) ==="
  vacuity_json="${report_dir}/j6_vacuity_probe.json"
  python3 "${script_dir}/j6_truth_vacuity_probe.py" \
    --bam "${P12_BAM}" \
    --truth-vcf "${P13_TRUTH_VCF}" \
    --regions-bed "${P13_REGIONS_BED}" \
    --eval-interval "${eval_interval}" \
    --min-depth "${J6_MIN_DEPTH:-5}" \
    --min-covered-truth "${J6_MIN_COVERED_TRUTH:-5}" \
    --json-out "${vacuity_json}"
else
  echo "=== stage NA12878_20k BAM ==="
  ./scripts/parity/realworld/02_stage_na12878_20k_bam.sh
fi

if [[ -n "${p12_json_backup}" ]]; then
  mv -f "${p12_json_backup}" "${p12_json}"
fi

if [[ "${J6_SKIP_HC:-0}" != "1" ]]; then
  echo ""
  echo "=== P12 scale: Java + Rust HC (${eval_interval}) ==="
  ./scripts/parity/run_p12_realworld_na12878_20k.sh
else
  echo ""
  echo "=== P12 scale skipped (J6_SKIP_HC=1) ==="
fi

echo ""
echo "=== P13 truth eval (stratified SNP/INDEL + gate) ==="
if [[ "${dense}" == "1" ]]; then
  p13_md_out="${report_dir}/p13_dense_truth_eval.md"
else
  p13_md_out="${repo_root}/parity/reports/p13_realworld_truth_eval.md"
fi
p13_args=(
  --java-vcf "${P13_JAVA_VCF:-${repo_root}/parity/reports/p12_realworld_na12878_20k.java.vcf}"
  --rust-vcf "${P13_RUST_VCF:-${repo_root}/parity/reports/p12_realworld_na12878_20k.rust.vcf}"
  --truth-vcf "${P13_TRUTH_VCF}"
  --json-out "${p13_json}"
  --md-out "${p13_md_out}"
  --chrom-filter "${P13_CHROM:-}"
  --eval-interval "${P13_EVAL_INTERVAL}"
  --thresholds-json "${thresholds}"
)
if [[ -f "${P13_REGIONS_BED}" ]]; then
  p13_args+=(--regions-bed "${P13_REGIONS_BED}")
fi
if [[ "${strict}" == "1" ]]; then
  p13_args+=(--strict-gate)
fi
python3 "${script_dir}/p13_truth_eval.py" "${p13_args[@]}"

j6_json="${report_dir}/j6_truth_summary.json"
j6_md="${report_dir}/j6_truth_summary.md"

echo ""
echo "=== L6 summary ==="
set +e
if [[ "${dense}" == "1" ]]; then
  parity_interval="${J6_PARITY_INTERVAL:-${eval_interval}}"
else
  parity_interval="${J6_PARITY_INTERVAL:-2:92300000-92350000}"
fi
python3 "${script_dir}/j6_truth_summarize.py" \
  --p12-json "${p12_json}" \
  --p13-json "${p13_json}" \
  --json-out "${j6_json}" \
  --md-out "${j6_md}" \
  --interval "${eval_interval}" \
  --parity-interval "${parity_interval}" \
  --thresholds "${thresholds}"
summary_exit=$?
set -e

cp -f "${log}" "${report_dir}/j6_truth_canonical.log"
cp -f "${j6_json}" "${report_dir}/j6_truth_canonical.json"
cp -f "${j6_md}" "${report_dir}/j6_truth_canonical.md"

if [[ "${summary_exit}" -ne 0 && "${strict}" == "1" ]]; then
  echo "[j6-truth] FAIL — see ${j6_md}" >&2
  exit "${summary_exit}"
fi

echo ""
echo "=== L6 hc-full-parity j6 truth PASS ==="
echo "summary: ${j6_md}"
echo "log: ${log}"
