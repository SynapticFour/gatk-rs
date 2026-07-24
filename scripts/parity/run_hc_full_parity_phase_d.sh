#!/usr/bin/env bash
# Phase D (D.1–D.4) — read path L1 gates. Runs Phase C (and thus B) first by default.
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

if [[ "${PARITY_SKIP_HC_FULL_PHASE_D:-0}" == "1" ]]; then
  echo "[hc-full-parity-phase-d] skipped (PARITY_SKIP_HC_FULL_PHASE_D=1)"
  exit 0
fi

if [[ "${PARITY_PHASE_D_SKIP_PHASE_C_CHECK:-0}" != "1" ]]; then
  echo "[hc-full-parity-phase-d] verifying Phase C first"
  PARITY_SKIP_HC_FULL_PHASE_C=0 "${repo_root}/scripts/parity/run_hc_full_parity_phase_c.sh"
fi

for gate in \
  run_hc_full_parity_d1_read_filters.sh \
  run_hc_full_parity_d2_downsample.sh \
  run_hc_full_parity_d2c_allele_biased.sh \
  run_hc_full_parity_d2c_contam.sh \
  run_hc_full_parity_d3_soft_clip.sh \
  run_hc_full_parity_d4_read_transform.sh; do
  echo "[hc-full-parity-phase-d] ${gate}"
  "${repo_root}/scripts/parity/${gate}"
done

echo "[hc-full-parity-phase-d] Phase D L1 gates (D.1–D.4): OK"
