#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
report_dir="${repo_root}/parity/reports"
mkdir -p "${report_dir}"
tmp_dir="${report_dir}/p3-region-diff-tmp"
mkdir -p "${tmp_dir}"

run_java="${repo_root}/scripts/parity/run_java_gatk.sh"
run_rust="${repo_root}/scripts/parity/run_rust_gatk.sh"

source_sam="${repo_root}/parity/fixtures/p3_optional_tags.sam"
source_bam="${tmp_dir}/p3_optional_tags.bam"
samtools view -bS "${source_sam}" > "${source_bam}"
samtools index "${source_bam}"

passed=0
failed=0
skipped=0
checks_json=()

extract_java_count() {
  local path="$1"
  python3 - <<PY
import re, pathlib
text = pathlib.Path("${path}").read_text(encoding="utf-8", errors="replace")
m = re.search(r"CountReads counted\\s+(\\d+)\\s+total reads", text)
if m:
    print(m.group(1))
else:
    m2 = re.search(r"Tool returned:\\s*(\\d+)", text, flags=re.MULTILINE)
    print(m2.group(1) if m2 else "")
PY
}

extract_rust_count() {
  local path="$1"
  python3 - <<PY
import re, pathlib
text = pathlib.Path("${path}").read_text(encoding="utf-8", errors="replace")
m = re.search(r"COUNT\\s*:\\s*(\\d+)", text)
print(m.group(1) if m else "")
PY
}

run_case() {
  local label="$1"
  local region="$2"
  local java_stdout="${report_dir}/${label}.java.stdout"
  local rust_stdout="${report_dir}/${label}.rust.stdout"
  local check_json="${report_dir}/${label}.json"

  set +e
  "${run_java}" "${java_stdout}" CountReads -I "${source_bam}" -L "${region}"
  java_exit=$?
  "${run_rust}" "${rust_stdout}" CountReadsInRegion -I "${source_bam}" -L "${region}"
  rust_exit=$?
  set -e

  if [[ "${java_exit}" -eq 127 && "${PARITY_ALLOW_MISSING_JAVA:-0}" == "1" ]]; then
    cat > "${check_json}" <<EOF
{
  "label": "${label}",
  "mode": "phase3-region-query-diff",
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

  if [[ "${java_exit}" -ne 0 || "${rust_exit}" -ne 0 ]]; then
    cat > "${check_json}" <<EOF
{
  "label": "${label}",
  "mode": "phase3-region-query-diff",
  "equal": false,
  "reason": "tool_exit_nonzero",
  "java_exit": ${java_exit},
  "rust_exit": ${rust_exit},
  "java_output": "${java_stdout}",
  "rust_output": "${rust_stdout}"
}
EOF
    failed=$((failed + 1))
    checks_json+=("${check_json}")
    return
  fi

  java_count="$(extract_java_count "${java_stdout}")"
  rust_count="$(extract_rust_count "${rust_stdout}")"
  equal=false
  if [[ -n "${java_count}" && -n "${rust_count}" && "${java_count}" == "${rust_count}" ]]; then
    equal=true
  fi

  cat > "${check_json}" <<EOF
{
  "label": "${label}",
  "mode": "phase3-region-query-diff",
  "equal": ${equal},
  "java_count": "${java_count}",
  "rust_count": "${rust_count}",
  "region": "${region}",
  "java_output": "${java_stdout}",
  "rust_output": "${rust_stdout}"
}
EOF

  if [[ "${equal}" == "true" ]]; then
    passed=$((passed + 1))
  else
    failed=$((failed + 1))
  fi
  checks_json+=("${check_json}")
}

run_case "p3-region-query-chr1-1-16" "chr1:1-16"
run_case "p3-region-query-chr1-17-32" "chr1:17-32"

summary_json="${report_dir}/p3-region-query-diff.json"
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
  echo "P3 region query differential checks failed: ${failed}"
  exit 1
fi

echo "P3 region query differential checks passed: ${passed}, skipped: ${skipped}"
