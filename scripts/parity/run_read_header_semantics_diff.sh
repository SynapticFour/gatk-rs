#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
report_dir="${repo_root}/parity/reports"
mkdir -p "${report_dir}"

run_java="${repo_root}/scripts/parity/run_java_gatk.sh"
run_rust="${repo_root}/scripts/parity/run_rust_gatk.sh"

passed=0
failed=0
skipped=0
checks_json=()

run_case() {
  local label="$1"
  local input_path="$2"
  local expected_ok="$3"
  local java_stdout="${report_dir}/${label}.java.stdout"
  local rust_stdout="${report_dir}/${label}.rust.stdout"
  local check_json="${report_dir}/${label}.json"

  set +e
  "${run_java}" "${java_stdout}" ValidateSamFile -I "${input_path}" -MODE SUMMARY
  java_exit=$?
  "${run_rust}" "${rust_stdout}" Validate "${input_path}" -t SAM
  rust_exit=$?
  set -e

  if [[ "${java_exit}" -eq 127 && "${PARITY_ALLOW_MISSING_JAVA:-0}" == "1" ]]; then
    cat > "${check_json}" <<EOF
{
  "label": "${label}",
  "mode": "read-header-semantics-runtime-diff",
  "equal": null,
  "skipped": true,
  "reason": "java_gatk_missing",
  "java_output": "${java_stdout}",
  "rust_output": "${rust_stdout}"
}
EOF
    skipped=$((skipped + 1))
    checks_json+=("${check_json}")
    return
  fi

  local java_ok=0
  local rust_ok=0
  [[ "${java_exit}" -eq 0 ]] && java_ok=1
  [[ "${rust_exit}" -eq 0 ]] && rust_ok=1
  local equal=0
  if [[ "${java_ok}" -eq "${rust_ok}" && "${rust_ok}" -eq "${expected_ok}" ]]; then
    equal=1
  fi

  cat > "${check_json}" <<EOF
{
  "label": "${label}",
  "mode": "read-header-semantics-runtime-diff",
  "equal": $( [[ "${equal}" -eq 1 ]] && echo "true" || echo "false" ),
  "expected_ok": ${expected_ok},
  "java_exit": ${java_exit},
  "rust_exit": ${rust_exit},
  "java_ok": $( [[ "${java_ok}" -eq 1 ]] && echo "true" || echo "false" ),
  "rust_ok": $( [[ "${rust_ok}" -eq 1 ]] && echo "true" || echo "false" ),
  "java_output": "${java_stdout}",
  "rust_output": "${rust_stdout}"
}
EOF

  if [[ "${equal}" -eq 1 ]]; then
    passed=$((passed + 1))
  else
    failed=$((failed + 1))
  fi
  checks_json+=("${check_json}")
}

run_case "read-header-semantics-valid-runtime-diff" \
  "${repo_root}/parity/fixtures/read_header_semantics_valid.sam" \
  1
run_case "read-header-semantics-missing-rg-runtime-diff" \
  "${repo_root}/parity/fixtures/read_header_semantics_missing_rg.sam" \
  0

summary_json="${report_dir}/read-header-semantics-diff.json"
{
  echo "{"
  echo "  \"passed\": ${passed},"
  echo "  \"failed\": ${failed},"
  echo "  \"skipped\": ${skipped},"
  echo "  \"checks\": ["
  for i in "${!checks_json[@]}"; do
    comma=","
    if [[ "${i}" -eq $((${#checks_json[@]} - 1)) ]]; then
      comma=""
    fi
    python3 - <<PY
import json, pathlib
print(json.dumps(json.loads(pathlib.Path("${checks_json[$i]}").read_text()), indent=2))
PY
    echo "${comma}"
  done
  echo "  ]"
  echo "}"
} > "${summary_json}"

if [[ "${failed}" -gt 0 ]]; then
  echo "Read header semantics runtime differential checks failed: ${failed}"
  exit 1
fi

echo "Read header semantics runtime differential checks passed: ${passed}, skipped: ${skipped}"
