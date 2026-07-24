#!/usr/bin/env bash
# L7-A3: second non-chr2 dense GIAB truth slice (chr21).
#
# Reuses J6_DENSE harness with a separate BAM OUT_DIR and report dir so chr20
# artifacts are not clobbered.
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"
cd "${repo_root}"

export J6_DENSE=1
export J6_DENSE_INTERVAL="${J6_DENSE_INTERVAL:-21:41200001-41250000}"
export J6_INTERVAL="${J6_INTERVAL:-${J6_DENSE_INTERVAL}}"
export J6_PARITY_INTERVAL="${J6_PARITY_INTERVAL:-${J6_INTERVAL}}"
export J6_DENSE_OUT_DIR="${J6_DENSE_OUT_DIR:-${repo_root}/parity/realworld/na12878_giab_window_chr21_b37}"
export J6_REPORT_DIR="${J6_REPORT_DIR:-${repo_root}/parity/reports/hc-full-parity-j6-dense-chr21}"
export J6_THRESHOLDS="${J6_THRESHOLDS:-${repo_root}/parity/fixtures/hc-full-parity/j6/thresholds_dense_chr21.json}"
export P12_REFERENCE="${P12_REFERENCE:-${repo_root}/parity/realworld/assets/hs37d5.simple.fa}"

# Reuse BAI already staged with the chr20 dense window when present.
chr20_bai="${repo_root}/parity/realworld/na12878_giab_window_b37/RMNISTHS_30xdownsample.bam.bai"
if [[ -z "${J6_DENSE_BAI_URL:-}" && -f "${chr20_bai}" ]]; then
  mkdir -p "${J6_DENSE_OUT_DIR}"
  if [[ ! -f "${J6_DENSE_OUT_DIR}/RMNISTHS_30xdownsample.bam.bai" ]]; then
    cp -f "${chr20_bai}" "${J6_DENSE_OUT_DIR}/RMNISTHS_30xdownsample.bam.bai"
  fi
fi

echo "=== L7 dense chr21 slice ==="
echo "interval=${J6_INTERVAL}"
echo "bam_out=${J6_DENSE_OUT_DIR}"
echo "reports=${J6_REPORT_DIR}"

exec ./scripts/parity/run_hc_full_parity_j6_truth.sh
