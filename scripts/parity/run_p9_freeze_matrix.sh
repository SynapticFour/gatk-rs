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
  echo "[p9-freeze] ${label}"
  "$@"
}

run_step "phase9-cli-contracts" ./scripts/parity/run_p9_cli_contracts.sh
run_step "phase9-hc-scaffold-diff" ./scripts/parity/run_p9_hc_scaffold_diff.sh
run_step "phase9-java-hc-smoke" ./scripts/parity/run_p9_java_hc_smoke.sh
run_step "phase9-mismatch-triage-check" ./scripts/parity/run_p9_mismatch_triage_check.sh
run_step "phase9-equivalence-summary" ./scripts/parity/build_p9_equivalence_summary.py

echo "P9 freeze matrix completed successfully."
