#!/usr/bin/env bash
# Phase H — GVCF / reference confidence L1 gates (+ P8 block contract).
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

if [[ "${PARITY_SKIP_HC_FULL_PHASE_H:-0}" == "1" ]]; then
  echo "[hc-full-parity-phase-h] skipped (PARITY_SKIP_HC_FULL_PHASE_H=1)"
  exit 0
fi

echo "[hc-full-parity-phase-h] run_hc_full_parity_h1_ref_confidence.sh"
"${repo_root}/scripts/parity/run_hc_full_parity_h1_ref_confidence.sh"

echo "[hc-full-parity-phase-h] run_hc_full_parity_h1_inactive.sh"
"${repo_root}/scripts/parity/run_hc_full_parity_h1_inactive.sh"

echo "[hc-full-parity-phase-h] run_hc_full_parity_h2_gvcf.sh"
"${repo_root}/scripts/parity/run_hc_full_parity_h2_gvcf.sh"

echo "[hc-full-parity-phase-h] run_p8_gvcf_contracts.sh"
"${repo_root}/scripts/parity/run_p8_gvcf_contracts.sh"

echo "[hc-full-parity-phase-h] p8_gvcf_block_diff_test (H.2.2 / L5 block semantics)"
cargo test -p gatk-haplotypecaller --test p8_gvcf_block_diff_test --locked

if [[ "${PARITY_SKIP_HC_FULL_DEFERRED_H:-0}" != "1" ]]; then
  echo "[hc-full-parity-phase-h] deferred H L5 scaffold"
  "${repo_root}/scripts/parity/run_hc_full_parity_h2_l5.sh"
fi

echo "[hc-full-parity-phase-h] Phase H L1 gates: OK"
