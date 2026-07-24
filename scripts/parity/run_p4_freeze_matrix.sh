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
  echo "[p4-freeze] ${label}"
  "$@"
}

# Step-62 freeze matrix: smoke profiles + Phase-4 activity surface + Java assembly-region harness + bench smoke.
run_step "smoke-profile-smoke" env PARITY_SMOKE_PROFILE=smoke ./scripts/parity/run_parity_smoke.sh
run_step "smoke-profile-extended" env PARITY_SMOKE_PROFILE=extended ./scripts/parity/run_parity_smoke.sh

run_step "phase2-haplotypecaller-lib" cargo test -p gatk-haplotypecaller --lib --locked
run_step "phase4-activity-region-creation-contract" cargo test -p gatk-haplotypecaller --test activity_region_creation_contract_test --locked
run_step "phase4-activity-profile-property-tests" cargo test -p gatk-haplotypecaller --test activity_profile_property_tests --locked
run_step "phase4-activity-repro-concurrency-contract" cargo test -p gatk-haplotypecaller --test p4_activity_repro_contract_test --locked
run_step "phase4-hc-assembly-region-interval-diff" ./scripts/parity/run_p4_active_region_interval_diff.sh
run_step "phase4-activity-profile-bench-smoke" cargo bench -p gatk-haplotypecaller --bench activity_profile --locked -- --quick

echo "P4 freeze matrix completed successfully."
