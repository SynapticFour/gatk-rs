#!/usr/bin/env bash
# Phase I — annotation manifest + core INFO plugin graph (parity v1).
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

if [[ "${PARITY_SKIP_HC_FULL_PHASE_I:-0}" == "1" ]]; then
  echo "[hc-full-parity-phase-i] skipped (PARITY_SKIP_HC_FULL_PHASE_I=1)"
  exit 0
fi

echo "[hc-full-parity-phase-i] run_hc_full_parity_i1_manifest.sh"
"${repo_root}/scripts/parity/run_hc_full_parity_i1_manifest.sh"

echo "[hc-full-parity-phase-i] run_hc_full_parity_i1_core.sh"
"${repo_root}/scripts/parity/run_hc_full_parity_i1_core.sh"

if [[ "${PARITY_SKIP_HC_FULL_DEFERRED_I:-0}" != "1" ]]; then
  echo "[hc-full-parity-phase-i] deferred I gates"
  "${repo_root}/scripts/parity/run_hc_full_parity_i1_standard.sh"
  "${repo_root}/scripts/parity/run_hc_full_parity_i1_as.sh"
  "${repo_root}/scripts/parity/run_hc_full_parity_i1_excess_het.sh"
  "${repo_root}/scripts/parity/run_hc_full_parity_i1_depth_hc.sh"
  "${repo_root}/scripts/parity/run_hc_full_parity_i1_plugins.sh"
fi

echo "[hc-full-parity-phase-i] Phase I L1 gates: OK"
