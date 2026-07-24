#!/usr/bin/env bash
# Runs 01 + 08–11, records exit codes, generates realworld_parity_evidence.{md,json}
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "${repo_root}"
rd="${repo_root}/scripts/parity/realworld"
manifest="${repo_root}/parity/reports/realworld_run_manifest.json"

run_one() {
  local phase="$1"
  local script="$2"
  local ec=0
  set +e
  "${rd}/${script}"
  ec=$?
  set -e
  python3 - <<PY
import json, pathlib
m = pathlib.Path("${manifest}")
rows = []
if m.exists():
    rows = json.loads(m.read_text(encoding="utf-8"))
rows.append({"phase": "${phase}", "script": "${script}", "exit_code": ${ec}})
m.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
PY
  if [[ "${ec}" -ne 0 ]]; then
    echo "run_foundation_evidence: ${script} FAILED exit=${ec}" >&2
    return "${ec}"
  fi
  return 0
}

echo "[]" > "${manifest}"

fail=0
run_one "01_env" "01_check_environment.sh" || fail=1
run_one "08_core_lib" "08_test_gatk_core_lib.sh" || fail=1
run_one "09_p3_io" "09_test_p3_io_conformance.sh" || fail=1
run_one "10_hc_lib" "10_test_haplotypecaller_lib.sh" || fail=1
run_one "11_smoke" "11_parity_smoke.sh" || fail=1

python3 "${rd}/generate_evidence_report.py" "${manifest}"

if [[ "${fail}" -ne 0 ]]; then
  echo "run_foundation_evidence: one or more steps failed" >&2
  exit 1
fi
echo "run_foundation_evidence: all steps OK"
exit 0
