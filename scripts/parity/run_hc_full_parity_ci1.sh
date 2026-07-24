#!/usr/bin/env bash
# CI.1 — hc-full-parity phase bundles B–J + strict L2 (220-case gate).
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"
export PARITY_HC_FULL_L2_STRICT="${PARITY_HC_FULL_L2_STRICT:-1}"
export PARITY_REQUIRE_SAMTOOLS="${PARITY_REQUIRE_SAMTOOLS:-1}"

if [[ "${PARITY_SKIP_HC_FULL_CI1:-0}" == "1" ]]; then
  echo "[hc-full-parity-ci1] skipped (PARITY_SKIP_HC_FULL_CI1=1)"
  exit 0
fi

for phase in b c d e f g h i j; do
  script="${repo_root}/scripts/parity/run_hc_full_parity_phase_${phase}.sh"
  echo "[hc-full-parity-ci1] ${script}"
  "${script}"
done

echo "[hc-full-parity-ci1] run_hc_full_parity_l2.sh (strict=${PARITY_HC_FULL_L2_STRICT})"
"${repo_root}/scripts/parity/run_hc_full_parity_l2.sh"

echo "[hc-full-parity-ci1] OK"
