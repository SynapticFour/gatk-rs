#!/usr/bin/env bash
# E-D02 — cyclic BAM gate (reuses e5 repeat fixtures as cycle proxy until dedicated BAM exists).
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

if [[ "${PARITY_SKIP_HC_FULL_E5_CYCLE:-0}" == "1" ]]; then
  echo "[hc-full-parity-e5-cycle] skipped (PARITY_SKIP_HC_FULL_E5_CYCLE=1)"
  exit 0
fi

echo "[hc-full-parity-e5-cycle] run_hc_full_parity_e5_assembly_graph_non_unique_summary.sh (cycle proxy)"
"${repo_root}/scripts/parity/run_hc_full_parity_e5_assembly_graph_non_unique_summary.sh"
