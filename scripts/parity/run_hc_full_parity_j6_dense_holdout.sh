#!/usr/bin/env bash
# L8 generalize: third dense GIAB holdout (chr20 offset window).
#
# Same contig as primary dense, different locus — required before deleting P12 pins
# (see docs/CLAIM_MATRIX.md).
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"
cd "${repo_root}"

export J6_DENSE=1
export J6_DENSE_INTERVAL="${J6_DENSE_INTERVAL:-20:15000000-15050000}"
export J6_INTERVAL="${J6_INTERVAL:-${J6_DENSE_INTERVAL}}"
export J6_PARITY_INTERVAL="${J6_PARITY_INTERVAL:-${J6_INTERVAL}}"
export J6_DENSE_OUT_DIR="${J6_DENSE_OUT_DIR:-${repo_root}/parity/realworld/na12878_giab_window_chr20_holdout_b37}"
export J6_REPORT_DIR="${J6_REPORT_DIR:-${repo_root}/parity/reports/hc-full-parity-j6-dense-holdout}"
export J6_THRESHOLDS="${J6_THRESHOLDS:-${repo_root}/parity/fixtures/hc-full-parity/j6/thresholds_dense_holdout.json}"
export P12_REFERENCE="${P12_REFERENCE:-${repo_root}/parity/realworld/assets/hs37d5.simple.fa}"

chr20_bai="${repo_root}/parity/realworld/na12878_giab_window_b37/RMNISTHS_30xdownsample.bam.bai"
if [[ -z "${J6_DENSE_BAI_URL:-}" && -f "${chr20_bai}" ]]; then
  mkdir -p "${J6_DENSE_OUT_DIR}"
  if [[ ! -f "${J6_DENSE_OUT_DIR}/RMNISTHS_30xdownsample.bam.bai" ]]; then
    cp -f "${chr20_bai}" "${J6_DENSE_OUT_DIR}/RMNISTHS_30xdownsample.bam.bai"
  fi
fi

echo "=== L8 dense chr20 holdout slice ==="
echo "interval=${J6_INTERVAL}"
echo "bam_out=${J6_DENSE_OUT_DIR}"
echo "reports=${J6_REPORT_DIR}"

exec ./scripts/parity/run_hc_full_parity_j6_truth.sh
