#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
report_dir="${repo_root}/parity/reports"
mkdir -p "${report_dir}"
cd "${repo_root}"

export LC_ALL=C
export TZ=UTC
export RAYON_NUM_THREADS="${RAYON_NUM_THREADS:-1}"
export PARITY_RANDOM_SEED="${PARITY_RANDOM_SEED:-1337}"
export PYTHONHASHSEED="${PYTHONHASHSEED:-${PARITY_RANDOM_SEED}}"
export SOURCE_DATE_EPOCH="${SOURCE_DATE_EPOCH:-1700000000}"

run_step() {
  local label="$1"
  shift
  echo "[p7-freeze] ${label}"
  "$@"
}

run_step "phase7-genotyping-contracts" ./scripts/parity/run_p7_genotyping_contracts.sh
run_step "phase7-genotype-field-diff" \
  cargo test -p gatk-haplotypecaller --test p7_genotype_field_diff_test --locked
run_step "phase7-mismatch-triage-check" ./scripts/parity/run_p7_mismatch_triage_check.sh
run_step "phase7-equivalence-summary" ./scripts/parity/build_p7_equivalence_summary.py

echo "P7 freeze matrix completed successfully."
