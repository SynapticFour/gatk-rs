#!/usr/bin/env bash
# Phase G — genotyping L1 gates (G.1 aggregate, G.2 PL + region genotype).
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

if [[ "${PARITY_SKIP_HC_FULL_PHASE_G:-0}" == "1" ]]; then
  echo "[hc-full-parity-phase-g] skipped (PARITY_SKIP_HC_FULL_PHASE_G=1)"
  exit 0
fi

echo "[hc-full-parity-phase-g] run_hc_full_parity_g1_genotyping_aggregate.sh"
"${repo_root}/scripts/parity/run_hc_full_parity_g1_genotyping_aggregate.sh"

echo "[hc-full-parity-phase-g] run_hc_full_parity_g2_pl.sh"
"${repo_root}/scripts/parity/run_hc_full_parity_g2_pl.sh"

echo "[hc-full-parity-phase-g] run_hc_full_parity_g2_region.sh"
"${repo_root}/scripts/parity/run_hc_full_parity_g2_region.sh"

if [[ "${PARITY_SKIP_HC_FULL_DEFERRED_G:-0}" != "1" ]]; then
  echo "[hc-full-parity-phase-g] deferred G gates"
  "${repo_root}/scripts/parity/run_hc_full_parity_g2_af.sh"
  "${repo_root}/scripts/parity/run_hc_full_parity_g3.sh"
  "${repo_root}/scripts/parity/run_hc_full_parity_g4.sh"
  "${repo_root}/scripts/parity/run_hc_full_parity_g4_force.sh"
  "${repo_root}/scripts/parity/run_hc_full_parity_g_subset.sh"
fi

echo "[hc-full-parity-phase-g] run_p7_genotyping_contracts.sh (unit contracts)"
"${repo_root}/scripts/parity/run_p7_genotyping_contracts.sh"

echo "[hc-full-parity-phase-g] Phase G L1 gates: OK"
