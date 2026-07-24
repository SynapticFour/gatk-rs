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
  echo "[p8-freeze] ${label}"
  "$@"
}

run_step "phase8-gvcf-contracts" ./scripts/parity/run_p8_gvcf_contracts.sh
run_step "phase8-gvcf-block-diff" \
  cargo test -p gatk-haplotypecaller --test p8_gvcf_block_diff_test --locked
run_step "phase8-mismatch-triage-check" ./scripts/parity/run_p8_mismatch_triage_check.sh
run_step "phase8-equivalence-summary" ./scripts/parity/build_p8_equivalence_summary.py

echo "P8 freeze matrix completed successfully."
