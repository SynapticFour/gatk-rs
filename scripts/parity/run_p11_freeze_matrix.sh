#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"

run_step() {
  local label="$1"
  shift
  echo "[p11-freeze] ${label}"
  "$@"
}

run_step "phase11-hc-output-activation-contracts" ./scripts/parity/run_p11_hc_output_activation_contracts.sh
run_step "phase11-hc-output-field-diff-smoke" ./scripts/parity/run_p11_hc_output_field_diff_smoke.sh
run_step "phase11-hc-output-field-diff-corpus" ./scripts/parity/run_p11_hc_output_field_diff_corpus.sh
run_step "phase11-mismatch-triage-check" ./scripts/parity/run_p11_mismatch_triage_check.sh
run_step "phase11-equivalence-summary" ./scripts/parity/build_p11_equivalence_summary.py

echo "P11 freeze matrix completed successfully."
