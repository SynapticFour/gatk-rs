#!/usr/bin/env bash
# Phase J — end-to-end VCF emission (assembly-region path) + j2 parity gates.
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

if [[ "${PARITY_SKIP_HC_FULL_PHASE_J:-0}" == "1" ]]; then
  echo "[hc-full-parity-phase-j] skipped (PARITY_SKIP_HC_FULL_PHASE_J=1)"
  exit 0
fi

echo "[hc-full-parity-phase-j] run_hc_full_parity_j2_cli.sh"
"${repo_root}/scripts/parity/run_hc_full_parity_j2_cli.sh"

echo "[hc-full-parity-phase-j] run_hc_full_parity_j2_vcf.sh"
"${repo_root}/scripts/parity/run_hc_full_parity_j2_vcf.sh"

echo "[hc-full-parity-phase-j] run_hc_full_parity_j2_format.sh"
"${repo_root}/scripts/parity/run_hc_full_parity_j2_format.sh"

if [[ "${PARITY_SKIP_HC_FULL_DEFERRED_J:-0}" != "1" ]]; then
  echo "[hc-full-parity-phase-j] deferred J gates"
  "${repo_root}/scripts/parity/run_hc_full_parity_j_modes.sh"
  "${repo_root}/scripts/parity/run_hc_full_parity_j_bamout.sh"
  "${repo_root}/scripts/parity/run_hc_full_parity_j_dragen.sh"
  "${repo_root}/scripts/parity/run_hc_full_parity_j_realworld.sh"
fi

echo "[hc-full-parity-phase-j] Phase J L1 gates: OK"
