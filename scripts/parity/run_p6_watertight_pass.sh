#!/usr/bin/env bash
set -euo pipefail

profile="${1:-lite}"
repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"

run_step() {
  local label="$1"
  shift
  echo "[p6-watertight:${profile}] ${label}"
  "$@"
}

case "${profile}" in
  lite)
    run_step "pairhmm-contracts" ./scripts/parity/run_p6_pairhmm_contracts.sh
    run_step "scalar-vector-equivalence" cargo test -p gatk-haplotypecaller --test p6_scalar_vector_equivalence_test --locked
    run_step "likelihood-vector-diff" ./scripts/parity/run_p6_likelihood_vector_diff.sh
    run_step "failure-mode-contracts" cargo test -p gatk-haplotypecaller --test p6_failure_mode_contract_test --locked
    run_step "equivalence-summary" ./scripts/parity/build_p6_equivalence_summary.py "${repo_root}"
    ;;
  live)
    run_step "wave-c-gates" ./scripts/parity/run_p6_wave_c_gates.sh
    run_step "determinism-matrix" ./scripts/parity/run_p6_determinism_matrix.sh
    run_step "mismatch-triage-check" ./scripts/parity/run_p6_mismatch_triage_check.sh
    run_step "equivalence-summary" ./scripts/parity/build_p6_equivalence_summary.py "${repo_root}"
    ;;
  full)
    run_step "full-haplotypecaller-tests" cargo test -p gatk-haplotypecaller --locked
    run_step "phase6-freeze-matrix" ./scripts/parity/run_p6_freeze_matrix.sh
    ;;
  *)
    echo "usage: $0 [lite|live|full]" >&2
    exit 2
    ;;
esac

echo "P6 watertight pass (${profile}) completed."
