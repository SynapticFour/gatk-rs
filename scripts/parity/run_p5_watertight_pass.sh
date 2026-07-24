#!/usr/bin/env bash
set -euo pipefail

profile="${1:-lite}"
repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"

run_step() {
  local label="$1"
  shift
  echo "[p5-watertight:${profile}] ${label}"
  "$@"
}

case "${profile}" in
  lite)
    run_step "runtime-candidate-diff" ./scripts/parity/run_p5_runtime_candidate_diff.sh
    run_step "determinism-matrix" ./scripts/parity/run_p5_determinism_matrix.sh
    run_step "mismatch-triage-check" ./scripts/parity/run_p5_mismatch_triage_check.sh
    run_step "regression-suite" cargo test -p gatk-haplotypecaller --test p5_assembly_regression_test --locked
    run_step "equivalence-summary" ./scripts/parity/build_p5_equivalence_summary.py "${repo_root}"
    ;;
  live)
    run_step "live-java-rust-candidate-diff" ./scripts/parity/run_p5_live_java_rust_diff.sh
    run_step "runtime-candidate-diff" ./scripts/parity/run_p5_runtime_candidate_diff.sh
    run_step "determinism-matrix" ./scripts/parity/run_p5_determinism_matrix.sh
    run_step "mismatch-triage-check" ./scripts/parity/run_p5_mismatch_triage_check.sh
    run_step "equivalence-summary" ./scripts/parity/build_p5_equivalence_summary.py "${repo_root}"
    ;;
  full)
    run_step "full-haplotypecaller-tests" cargo test -p gatk-haplotypecaller --locked
    run_step "phase5-freeze-matrix" ./scripts/parity/run_p5_freeze_matrix.sh
    ;;
  *)
    echo "usage: $0 [lite|live|full]" >&2
    exit 2
    ;;
esac

echo "P5 watertight pass (${profile}) completed."
