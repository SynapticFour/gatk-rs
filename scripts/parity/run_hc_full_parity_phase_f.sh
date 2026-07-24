#!/usr/bin/env bash
# Phase F (F.1+) — PairHMM / ReadLikelihoods L1 gates.
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

if [[ "${PARITY_SKIP_HC_FULL_PHASE_F:-0}" == "1" ]]; then
  echo "[hc-full-parity-phase-f] skipped (PARITY_SKIP_HC_FULL_PHASE_F=1)"
  exit 0
fi

echo "[hc-full-parity-phase-f] run_hc_full_parity_f1_pairhmm_likelihoods.sh"
"${repo_root}/scripts/parity/run_hc_full_parity_f1_pairhmm_likelihoods.sh"

echo "[hc-full-parity-phase-f] run_hc_full_parity_f2_pairhmm_native.sh"
"${repo_root}/scripts/parity/run_hc_full_parity_f2_pairhmm_native.sh"

echo "[hc-full-parity-phase-f] run_hc_full_parity_f3 (if present)"
if [[ -f "${repo_root}/scripts/parity/run_hc_full_parity_f3.sh" ]]; then
  "${repo_root}/scripts/parity/run_hc_full_parity_f3.sh" || true
fi

echo "[hc-full-parity-phase-f] Phase F L1 gates (F.1 + F.2): OK"
