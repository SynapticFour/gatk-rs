#!/usr/bin/env bash
# Phase 3 unmapped + supplementary matrix: Java CountReads vs Rust CountReadsInRegion.
# Fixture: primary r1 chr1:100-129 (30M), supplementary r1 chr1:200-209 (10M), fully unmapped r_unmap.
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
report_dir="${repo_root}/parity/reports"
mkdir -p "${report_dir}"
tmp_dir="${report_dir}/p3-unmapped-supp-diff-tmp"
mkdir -p "${tmp_dir}"

run_java="${repo_root}/scripts/parity/run_java_gatk.sh"
run_rust="${repo_root}/scripts/parity/run_rust_gatk.sh"

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
  local source_bam="$2"
  local fixture_name="$3"
  local region="$4"
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
  "mode": "phase3-unmapped-supplementary-runtime-diff",
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
  "mode": "phase3-unmapped-supplementary-runtime-diff",
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
  "mode": "phase3-unmapped-supplementary-runtime-diff",
  "equal": ${equal},
  "java_count": "${java_count}",
  "rust_count": "${rust_count}",
  "fixture": "${fixture_name}",
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

prepare_fixture_bam() {
  local fixture_file="$1"
  local out_bam="$2"
  samtools view -bS "${fixture_file}" > "${out_bam}"
  samtools index "${out_bam}"
}

run_fixture_cases() {
  local fixture_name="$1"
  local primary_window="$2"
  local supplementary_window="$3"
  local union_window="$4"
  local empty_window="$5"
  local fixture_file="${repo_root}/parity/fixtures/${fixture_name}"
  local source_bam="${tmp_dir}/${fixture_name%.sam}.bam"
  prepare_fixture_bam "${fixture_file}" "${source_bam}"

  run_case "p3-unmapped-supp-${fixture_name%.sam}-primary-only" "${source_bam}" "${fixture_name}" "${primary_window}"
  run_case "p3-unmapped-supp-${fixture_name%.sam}-supp-only" "${source_bam}" "${fixture_name}" "${supplementary_window}"
  run_case "p3-unmapped-supp-${fixture_name%.sam}-both-segments" "${source_bam}" "${fixture_name}" "${union_window}"
  run_case "p3-unmapped-supp-${fixture_name%.sam}-empty-window" "${source_bam}" "${fixture_name}" "${empty_window}"
}

# Dataset 1: base fixture
run_fixture_cases "p3_unmapped_supplementary.sam" "chr1:95-125" "chr1:195-210" "chr1:95-210" "chr1:1-10"
# Dataset 2: alternative fixture with shifted windows
run_fixture_cases "p3_unmapped_supplementary_alt.sam" "chr1:390-420" "chr1:695-715" "chr1:390-715" "chr1:1-50"

summary_json="${report_dir}/p3-unmapped-supplementary-diff.json"
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
  echo "P3 unmapped/supplementary differential checks failed: ${failed}"
  exit 1
fi

echo "P3 unmapped/supplementary differential checks passed: ${passed}, skipped: ${skipped}"
