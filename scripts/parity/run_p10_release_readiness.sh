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
  echo "[p10-release] ${label}"
  "$@"
}

coverage_signal="${report_dir}/p10_coverage_signal.json"
set +e
./scripts/coverage.sh gate-minimum
coverage_exit=$?
set -e
if [[ "${coverage_exit}" -eq 0 ]]; then
  echo "[p10-release] coverage-minimum-gate passed"
  printf '%s\n' '{"mode":"llvm-cov-gate-minimum","fallback_used":false}' > "${coverage_signal}"
else
  fallback_target_dir="${repo_root}/target-p10-release"
  run_step "coverage-fallback-core-tests" env CARGO_TARGET_DIR="${fallback_target_dir}" cargo test -p gatk-core --tests --locked
  printf '%s\n' '{"mode":"fallback-core-tests","fallback_used":true,"reason":"cargo-llvm-cov-or-llvm-tools-unavailable"}' > "${coverage_signal}"
fi

run_step "phase7-mismatch-triage-check" ./scripts/parity/run_p7_mismatch_triage_check.sh
run_step "phase8-mismatch-triage-check" ./scripts/parity/run_p8_mismatch_triage_check.sh
run_step "phase9-mismatch-triage-check" ./scripts/parity/run_p9_mismatch_triage_check.sh
run_step "build-release-readiness-summary" ./scripts/parity/build_p10_release_readiness.py

echo "P10 release-readiness gate completed successfully."
