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
  local expected_ok="$2"
  local java_args="$3"
  local rust_args="$4"

  local java_stdout="${report_dir}/${label}.java.stdout"
  local rust_stdout="${report_dir}/${label}.rust.stdout"
  local check_json="${report_dir}/${label}.json"

  set +e
  # shellcheck disable=SC2206
  local java_arr=( ${java_args} )
  # shellcheck disable=SC2206
  local rust_arr=( ${rust_args} )
  "${run_java}" "${java_stdout}" "${java_arr[@]}"
  java_exit=$?
  "${run_rust}" "${rust_stdout}" "${rust_arr[@]}"
  rust_exit=$?
  set -e

  if [[ "${java_exit}" -eq 127 && "${PARITY_ALLOW_MISSING_JAVA:-0}" == "1" ]]; then
    cat > "${check_json}" <<EOF
{
  "label": "${label}",
  "mode": "phase3-truncation-corruption-runtime-diff",
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

  local java_ok=false
  local rust_ok=false
  [[ "${java_exit}" -eq 0 ]] && java_ok=true
  [[ "${rust_exit}" -eq 0 ]] && rust_ok=true

  local equal=false
  if [[ "${java_ok}" == "${rust_ok}" ]]; then
    if [[ "${expected_ok}" -eq 1 && "${java_ok}" == "true" ]]; then
      equal=true
    elif [[ "${expected_ok}" -eq 0 && "${java_ok}" == "false" ]]; then
      equal=true
    fi
  fi

  cat > "${check_json}" <<EOF
{
  "label": "${label}",
  "mode": "phase3-truncation-corruption-runtime-diff",
  "equal": ${equal},
  "expected_ok": ${expected_ok},
  "java_exit": ${java_exit},
  "rust_exit": ${rust_exit},
  "java_ok": ${java_ok},
  "rust_ok": ${rust_ok},
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

reference="${repo_root}/parity/fixtures/reference.fa"

# Valid controls.
run_case "p3-trunc-corr-control-sam-valid" 1 \
  "ValidateSamFile -I ${repo_root}/parity/fixtures/sample.sam -MODE SUMMARY -IGNORE_WARNINGS true -REFERENCE_SEQUENCE ${reference}" \
  "Validate ${repo_root}/parity/fixtures/sample.sam -t SAM -R ${reference}"
run_case "p3-trunc-corr-control-vcf-valid" 1 \
  "ValidateVariants -V ${repo_root}/parity/fixtures/sample.vcf -R ${reference}" \
  "Validate ${repo_root}/parity/fixtures/sample.vcf -t VCF -R ${reference}"

# Truncation + corruption matrix.
run_case "p3-trunc-corr-bam-truncated-header" 0 \
  "ValidateSamFile -I ${repo_root}/parity/fixtures/p3_truncated_header.bam -MODE SUMMARY" \
  "Validate ${repo_root}/parity/fixtures/p3_truncated_header.bam -t BAM"
run_case "p3-trunc-corr-bam-bad-magic" 0 \
  "ValidateSamFile -I ${repo_root}/parity/fixtures/p3_malformed_bad_magic.bam -MODE SUMMARY" \
  "Validate ${repo_root}/parity/fixtures/p3_malformed_bad_magic.bam -t BAM"
run_case "p3-trunc-corr-sam-short-record" 0 \
  "ValidateSamFile -I ${repo_root}/parity/fixtures/p3_malformed_short_record.sam -MODE SUMMARY -REFERENCE_SEQUENCE ${reference}" \
  "Validate ${repo_root}/parity/fixtures/p3_malformed_short_record.sam -t SAM -R ${reference}"
run_case "p3-trunc-corr-sam-cigar-seq-mismatch" 0 \
  "ValidateSamFile -I ${repo_root}/parity/fixtures/p3_malformed_cigar_seq_mismatch.sam -MODE SUMMARY -REFERENCE_SEQUENCE ${reference}" \
  "Validate ${repo_root}/parity/fixtures/p3_malformed_cigar_seq_mismatch.sam -t SAM -R ${reference}"
run_case "p3-trunc-corr-vcf-short-record" 0 \
  "ValidateVariants -V ${repo_root}/parity/fixtures/p3_malformed_short_record.vcf -R ${reference}" \
  "Validate ${repo_root}/parity/fixtures/p3_malformed_short_record.vcf -t VCF -R ${reference}"
run_case "p3-trunc-corr-vcf-bad-pos" 0 \
  "ValidateVariants -V ${repo_root}/parity/fixtures/p3_malformed_bad_pos.vcf -R ${reference}" \
  "Validate ${repo_root}/parity/fixtures/p3_malformed_bad_pos.vcf -t VCF -R ${reference}"

summary_json="${report_dir}/p3-truncation-corruption-diff.json"
{
  echo "{"
  echo "  \"passed\": ${passed},"
  echo "  \"failed\": ${failed},"
  echo "  \"skipped\": ${skipped},"
  echo "  \"checks\": ["
  for i in "${!checks_json[@]}"; do
    comma=",";
    if [[ "${i}" -eq $((${#checks_json[@]} - 1)) ]]; then
      comma=""
    fi
    python3 - <<PY2
import json, pathlib
print(json.dumps(json.loads(pathlib.Path("${checks_json[$i]}").read_text()), indent=2))
PY2
    echo "${comma}"
  done
  echo "  ]"
  echo "}"
} > "${summary_json}"

if [[ "${failed}" -gt 0 ]]; then
  echo "P3 truncation/corruption differential checks failed: ${failed}"
  exit 1
fi

echo "P3 truncation/corruption differential checks passed: ${passed}, skipped: ${skipped}"
