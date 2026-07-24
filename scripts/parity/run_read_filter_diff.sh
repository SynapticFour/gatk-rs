#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
report_dir="${repo_root}/parity/reports"
mkdir -p "${report_dir}"
tmp_dir="${report_dir}/read-filter-diff-tmp"
mkdir -p "${tmp_dir}"

run_java="${repo_root}/scripts/parity/run_java_gatk.sh"
run_rust="${repo_root}/scripts/parity/run_rust_gatk.sh"
compare_sam_py="${repo_root}/scripts/parity/compare_sam_parity.py"

passed=0
failed=0
skipped=0
checks_json=()

run_case() {
  local label="$1"
  local input_path="$2"
  local java_stdout="${report_dir}/${label}.java.stdout"
  local rust_stdout="${report_dir}/${label}.rust.stdout"
  local java_sam="${report_dir}/${label}.java.sam"
  local rust_sam="${report_dir}/${label}.rust.sam"
  local check_json="${report_dir}/${label}.json"

  set +e
  "${run_java}" "${java_stdout}" PrintReads \
    -I "${input_path}" \
    -O "${java_sam}" \
    --read-filter MappedReadFilter \
    --read-filter NotSecondaryAlignmentReadFilter \
    --read-filter NotSupplementaryAlignmentReadFilter \
    --read-filter MappingQualityReadFilter \
    --minimum-mapping-quality 20 \
    --read-filter NotDuplicateReadFilter
  java_exit=$?
  "${run_rust}" "${rust_stdout}" FilterReads \
    -I "${input_path}" \
    -O "${rust_sam}" \
    --min-mapping-quality 20
  rust_exit=$?
  set -e

  if [[ "${java_exit}" -eq 127 && "${PARITY_ALLOW_MISSING_JAVA:-0}" == "1" ]]; then
    cat > "${check_json}" <<EOF
{
  "label": "${label}",
  "mode": "read-filter-runtime-diff",
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
  "mode": "read-filter-runtime-diff",
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

  set +e
  python3 "${compare_sam_py}" \
    --java-sam "${java_sam}" \
    --rust-sam "${rust_sam}" \
    --label "${label}" \
    --json-out "${check_json}"
  cmp_exit=$?
  set -e
  if [[ "${cmp_exit}" -eq 0 ]]; then
    passed=$((passed + 1))
  else
    failed=$((failed + 1))
  fi
  checks_json+=("${check_json}")
}

make_bam_slice() {
  local source_bam="$1"
  local start_record="$2"
  local record_count="$3"
  local out_bam="$4"
  local out_sam="${out_bam%.bam}.sam"
  local end_record=$((start_record + record_count))

  samtools view -H "${source_bam}" > "${out_sam}"
  samtools view "${source_bam}" | awk -v s="${start_record}" -v e="${end_record}" 'NR >= s && NR < e' >> "${out_sam}"
  samtools view -bS "${out_sam}" > "${out_bam}"
}

# Step 33 target: synthetic matrix fixture
run_case "read-filter-synthetic-runtime-diff" "${repo_root}/parity/fixtures/read_filter_slice.sam"

# Step 34 target: BAM fixture differential checks on multiple small slices.
source_bam="${tmp_dir}/read-filter-source.bam"
samtools view -bS "${repo_root}/parity/fixtures/read_filter_slice.sam" > "${source_bam}"
slice_a="${tmp_dir}/bam-slice-a.bam"
slice_b="${tmp_dir}/bam-slice-b.bam"
make_bam_slice "${source_bam}" 1 3 "${slice_a}"
make_bam_slice "${source_bam}" 4 3 "${slice_b}"
run_case "read-filter-bam-slice-a-runtime-diff" "${slice_a}"
run_case "read-filter-bam-slice-b-runtime-diff" "${slice_b}"

# Optional extra fixture from core test data when a usable BAM is present.
core_bam="${repo_root}/gatk-core/src/tests/test_data/sample.bam"
if [[ -f "${core_bam}" ]] && samtools view "${core_bam}" >/dev/null 2>&1; then
  core_slice="${tmp_dir}/core-bam-slice.bam"
  make_bam_slice "${core_bam}" 1 20 "${core_slice}"
  run_case "read-filter-core-bam-slice-runtime-diff" "${core_slice}"
fi

summary_json="${report_dir}/read-filter-diff.json"
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
  echo "Read-filter runtime differential checks failed: ${failed}"
  exit 1
fi

echo "Read-filter runtime differential checks passed: ${passed}, skipped: ${skipped}"
