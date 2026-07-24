#!/usr/bin/env bash
# Phase C (C.1–C.5) — activity / isActive + genotype-likelihood (MinimalGenotyping path) L1 gates.
# Requires Phase B (see docs/CLAIM_MATRIX.md).
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

if [[ "${PARITY_SKIP_HC_FULL_PHASE_C:-0}" == "1" ]]; then
  echo "[hc-full-parity-phase-c] skipped (PARITY_SKIP_HC_FULL_PHASE_C=1)"
  exit 0
fi

if [[ "${PARITY_PHASE_C_SKIP_PHASE_B_CHECK:-0}" != "1" ]]; then
  echo "[hc-full-parity-phase-c] verifying Phase B first"
  PARITY_SKIP_HC_FULL_PHASE_B=0 "${repo_root}/scripts/parity/run_hc_full_parity_phase_b.sh"
fi

for gate in \
  run_hc_full_parity_c1_raw_activity.sh \
  run_hc_full_parity_c2_smoothed_activity.sh \
  run_hc_full_parity_c3_active_locus.sh \
  run_hc_full_parity_c4_gl.sh \
  run_hc_full_parity_c5_multi.sh \
  run_hc_full_parity_c5_force.sh; do
  echo "[hc-full-parity-phase-c] ${gate}"
  "${repo_root}/scripts/parity/${gate}"
done

echo "[hc-full-parity-phase-c] Phase C L1 gates (C.1–C.5): OK"
