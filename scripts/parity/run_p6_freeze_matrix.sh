#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
report_dir="${repo_root}/parity/reports"
mkdir -p "${report_dir}"

export LC_ALL=C
export TZ=UTC
export RAYON_NUM_THREADS="${RAYON_NUM_THREADS:-1}"
export PARITY_RANDOM_SEED="${PARITY_RANDOM_SEED:-1337}"
export PYTHONHASHSEED="${PYTHONHASHSEED:-${PARITY_RANDOM_SEED}}"
export SOURCE_DATE_EPOCH="${SOURCE_DATE_EPOCH:-1700000000}"

run_step() {
  local label="$1"
  shift
  echo "[p6-freeze] ${label}"
  "$@"
}

run_step "phase6-pairhmm-contracts" ./scripts/parity/run_p6_pairhmm_contracts.sh
run_step "phase6-scalar-vector-equivalence" cargo test -p gatk-haplotypecaller --test p6_scalar_vector_equivalence_test --locked
run_step "phase6-likelihood-vector-diff" ./scripts/parity/run_p6_likelihood_vector_diff.sh
run_step "phase6-boundary-artifact-contracts" cargo test -p gatk-haplotypecaller --test p6_boundary_artifact_contract_test --locked
run_step "phase6-fp-policy-contract" cargo test -p gatk-haplotypecaller --test p6_fp_policy_contract_test --locked
run_step "phase6-failure-mode-contracts" cargo test -p gatk-haplotypecaller --test p6_failure_mode_contract_test --locked
run_step "phase6-determinism-matrix" ./scripts/parity/run_p6_determinism_matrix.sh
run_step "phase6-mismatch-triage-check" ./scripts/parity/run_p6_mismatch_triage_check.sh
run_step "phase6-pairhmm-bench-smoke" env CARGO_PROFILE_BENCH_OPT_LEVEL=1 cargo bench -p gatk-haplotypecaller --bench pairhmm --locked -- --quick
run_step "phase6-equivalence-summary" ./scripts/parity/build_p6_equivalence_summary.py

echo "P6 freeze matrix completed successfully."
