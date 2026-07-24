#!/usr/bin/env bash
# L11 production sign-off gate battery (extends L9: dense FORMAT+F1 + optional P12).
#
# Assumes release HC binary already built. Prefer:
#   source scripts/parity/m4_disk_guard.sh && m4_require_free_gb 12
#   export CARGO_TARGET_DIR="$PWD/target" CARGO_BUILD_JOBS=1 RAYON_NUM_THREADS=2
#   ./scripts/parity/build_gatk_rs_release.sh
#   L9_REGEN_RUST=1 L9_RUN_P12=1 ./scripts/parity/run_l11_signoff_gates.sh
#   ./scripts/parity/run_p12_l4_signoff.sh
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"
cd "${repo_root}"

export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$PWD/target}"
export CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-1}"
export RAYON_NUM_THREADS="${RAYON_NUM_THREADS:-2}"
export P12_SKIP_JAVA="${P12_SKIP_JAVA:-1}"

ts="$(date -u +%Y%m%dT%H%M%SZ)"
log="${repo_root}/parity/reports/l11_signoff_gates_${ts}.log"
mkdir -p "${repo_root}/parity/reports"

{
  echo "=== L11 sign-off gates ${ts} ==="
  echo "CARGO_TARGET_DIR=${CARGO_TARGET_DIR}"
  python3 "${repo_root}/scripts/parity/excellence_gates_audit.py"
  # Reuse L9 battery (chr20/chr21/holdout FORMAT+F1; optional L9_RUN_P12=1 → L3)
  ./scripts/parity/run_l9_signoff_gates.sh
  echo ""
  echo "=== L11 sign-off gates PASS ==="
  echo "log: ${log}"
  echo "Next: ./scripts/parity/run_p12_l4_signoff.sh"
  echo "See docs/CLAIM_MATRIX.md"
} 2>&1 | tee "${log}"
