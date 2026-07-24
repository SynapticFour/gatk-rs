#!/usr/bin/env bash
# Regenerate equivalence_report.* from an existing OUT_DIR and assert strict gates (no Docker).
# Usage:
#   ./scripts/parity/realworld/pipeline/run_realworld_equivalence_selfcheck.sh
#   OUT_DIR=/path/to/run ./scripts/parity/realworld/pipeline/run_realworld_equivalence_selfcheck.sh
set -euo pipefail

_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "${_SCRIPT_DIR}/common.sh"

OUT="${OUT_DIR}"
REPORT_PY="${_SCRIPT_DIR}/realworld_equivalence_report.py"
ASSERT_PY="${_SCRIPT_DIR}/assert_realworld_equivalence.py"

echo "[selfcheck] OUT_DIR=${OUT}"
python3 "${REPORT_PY}" "${OUT}"
python3 "${ASSERT_PY}" "${OUT}"
echo "[selfcheck] PASS — report regenerated and strict JSON gates satisfied."
