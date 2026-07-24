#!/usr/bin/env bash
# Phase PRE — read preparation L1 gates (after Phase D).
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

if [[ "${PARITY_SKIP_HC_FULL_PHASE_PRE:-0}" == "1" ]]; then
  echo "[hc-full-parity-phase-pre] skipped (PARITY_SKIP_HC_FULL_PHASE_PRE=1)"
  exit 0
fi

if [[ "${PARITY_PHASE_PRE_SKIP_PHASE_D_CHECK:-0}" != "1" ]]; then
  echo "[hc-full-parity-phase-pre] verifying Phase D first"
  PARITY_PHASE_D_SKIP_PHASE_C_CHECK=1 PARITY_SKIP_HC_FULL_PHASE_D=0 \
    "${repo_root}/scripts/parity/run_hc_full_parity_phase_d.sh"
fi

for gate in run_hc_full_parity_pre_unclip.sh run_hc_full_parity_pre_len.sh run_hc_full_parity_pre_mq.sh run_hc_full_parity_pre_overlap.sh; do
  echo "[hc-full-parity-phase-pre] ${gate}"
  "${repo_root}/scripts/parity/${gate}"
done

if [[ "${PARITY_SKIP_HC_FULL_DEFERRED_PRE:-0}" != "1" ]]; then
  echo "[hc-full-parity-phase-pre] deferred PRE-D01"
  "${repo_root}/scripts/parity/run_hc_full_parity_pre_dragstr.sh"
fi

echo "[hc-full-parity-phase-pre] Phase PRE L1 gates (PRE.1–PRE.4): OK"
