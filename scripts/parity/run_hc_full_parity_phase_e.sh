#!/usr/bin/env bash
# Phase E (E.1+) — assembly graph L1 gates.
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

if [[ "${PARITY_SKIP_HC_FULL_PHASE_E:-0}" == "1" ]]; then
  echo "[hc-full-parity-phase-e] skipped (PARITY_SKIP_HC_FULL_PHASE_E=1)"
  exit 0
fi

echo "[hc-full-parity-phase-e] run_hc_full_parity_e0_assemble.sh"
"${repo_root}/scripts/parity/run_hc_full_parity_e0_assemble.sh"

echo "[hc-full-parity-phase-e] run_hc_full_parity_e1_assembly_graph.sh"
"${repo_root}/scripts/parity/run_hc_full_parity_e1_assembly_graph.sh"

echo "[hc-full-parity-phase-e] run_hc_full_parity_e1_rec.sh"
"${repo_root}/scripts/parity/run_hc_full_parity_e1_rec.sh"

echo "[hc-full-parity-phase-e] run_hc_full_parity_e2_assembly_graph_multi.sh"
"${repo_root}/scripts/parity/run_hc_full_parity_e2_assembly_graph_multi.sh"

echo "[hc-full-parity-phase-e] run_hc_full_parity_e3_assembly_graph_summary.sh"
"${repo_root}/scripts/parity/run_hc_full_parity_e3_assembly_graph_summary.sh"

echo "[hc-full-parity-phase-e] run_hc_full_parity_e4_assembly_graph_dangling_summary.sh"
"${repo_root}/scripts/parity/run_hc_full_parity_e4_assembly_graph_dangling_summary.sh"

echo "[hc-full-parity-phase-e] run_hc_full_parity_e5_assembly_graph_non_unique_summary.sh"
"${repo_root}/scripts/parity/run_hc_full_parity_e5_assembly_graph_non_unique_summary.sh"

echo "[hc-full-parity-phase-e] run_hc_full_parity_e6_assembly_haplotype_cigars.sh"
"${repo_root}/scripts/parity/run_hc_full_parity_e6_assembly_haplotype_cigars.sh"

echo "[hc-full-parity-phase-e] run_hc_full_parity_e7_assembly_haplotypes.sh"
"${repo_root}/scripts/parity/run_hc_full_parity_e7_assembly_haplotypes.sh"

for gate in run_hc_full_parity_e7_kbest.sh run_hc_full_parity_e7_cap.sh \
  run_hc_full_parity_e7_artificial.sh run_hc_full_parity_e7_junction.sh \
  run_hc_full_parity_e7_edges.sh; do
  echo "[hc-full-parity-phase-e] ${gate}"
  "${repo_root}/scripts/parity/${gate}"
done

echo "[hc-full-parity-phase-e] run_hc_full_parity_e8_seqgraph.sh"
"${repo_root}/scripts/parity/run_hc_full_parity_e8_seqgraph.sh"

echo "[hc-full-parity-phase-e] run_hc_full_parity_e2e.sh"
"${repo_root}/scripts/parity/run_hc_full_parity_e2e.sh"

if [[ "${PARITY_SKIP_HC_FULL_DEFERRED_E:-0}" != "1" ]]; then
  echo "[hc-full-parity-phase-e] deferred E gates"
  "${repo_root}/scripts/parity/run_hc_full_parity_e_debug.sh"
  "${repo_root}/scripts/parity/run_hc_full_parity_e5_cycle.sh"
fi

echo "[hc-full-parity-phase-e] Phase E L1 gates (E.0–E.8 + E.7.x + E2E): OK"
