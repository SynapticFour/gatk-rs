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
  echo "[p5-freeze] ${label}"
  "$@"
}

# Step-76 freeze bundle for local assembly graph.
run_step "phase5-assembly-core-tests" cargo test -p gatk-haplotypecaller --lib assembly::tests:: --locked
run_step "phase5-runtime-candidate-diff" ./scripts/parity/run_p5_runtime_candidate_diff.sh
run_step "phase5-haplotype-candidate-diff" ./scripts/parity/run_p5_haplotype_candidate_diff.sh
run_step "phase5-assembly-regression-suite" cargo test -p gatk-haplotypecaller --test p5_assembly_regression_test --locked
run_step "phase5-determinism-matrix" ./scripts/parity/run_p5_determinism_matrix.sh
run_step "phase5-assembly-stability-contract" ./scripts/parity/run_p5_assembly_stability_contract.sh
run_step "phase5-mismatch-triage-check" ./scripts/parity/run_p5_mismatch_triage_check.sh
run_step "phase5-assembly-bench-smoke" env CARGO_PROFILE_BENCH_OPT_LEVEL=1 cargo bench -p gatk-haplotypecaller --bench assembly_graph --locked -- --quick
run_step "phase5-equivalence-summary" ./scripts/parity/build_p5_equivalence_summary.py

echo "P5 freeze matrix completed successfully."
