#!/usr/bin/env bash
# All deferred parity items (G/H/I/J/E/PRE) — L1 scaffold gates.
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

if [[ "${PARITY_SKIP_HC_FULL_DEFERRED:-0}" == "1" ]]; then
  echo "[hc-full-parity-deferred] skipped (PARITY_SKIP_HC_FULL_DEFERRED=1)"
  exit 0
fi

run() {
  echo "[hc-full-parity-deferred] $*"
  "$@"
}

# Phase G deferred
run "${repo_root}/scripts/parity/run_hc_full_parity_g2_af.sh"
run "${repo_root}/scripts/parity/run_hc_full_parity_g3.sh"
run "${repo_root}/scripts/parity/run_hc_full_parity_g4.sh"
run "${repo_root}/scripts/parity/run_hc_full_parity_g4_force.sh"
run "${repo_root}/scripts/parity/run_hc_full_parity_g_subset.sh"

# Phase H deferred (L5 scaffold; h2-blocks remains in phase_h)
run "${repo_root}/scripts/parity/run_hc_full_parity_h2_l5.sh"

# Phase I deferred
run "${repo_root}/scripts/parity/run_hc_full_parity_i1_standard.sh"
run "${repo_root}/scripts/parity/run_hc_full_parity_i1_as.sh"
run "${repo_root}/scripts/parity/run_hc_full_parity_i1_excess_het.sh"
run "${repo_root}/scripts/parity/run_hc_full_parity_i1_depth_hc.sh"
run "${repo_root}/scripts/parity/run_hc_full_parity_i1_plugins.sh"

# Phase J deferred
run "${repo_root}/scripts/parity/run_hc_full_parity_j_modes.sh"
run "${repo_root}/scripts/parity/run_hc_full_parity_j_bamout.sh"
run "${repo_root}/scripts/parity/run_hc_full_parity_j_dragen.sh"
run "${repo_root}/scripts/parity/run_hc_realworld_parity.sh"

# Phase E deferred
run "${repo_root}/scripts/parity/run_hc_full_parity_e_debug.sh"
run "${repo_root}/scripts/parity/run_hc_full_parity_e5_cycle.sh"

# PRE deferred
run "${repo_root}/scripts/parity/run_hc_full_parity_pre_dragstr.sh"

echo "[hc-full-parity-deferred] all deferred L1 gates: OK"
