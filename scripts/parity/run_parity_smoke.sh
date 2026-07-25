#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
report_dir="${repo_root}/parity/reports"
mkdir -p "${report_dir}"
run_tmp="${report_dir}/tmp-run"
mkdir -p "${run_tmp}"
export TMPDIR="${TMPDIR:-${run_tmp}}"
export RAYON_NUM_THREADS="${RAYON_NUM_THREADS:-1}"

export LC_ALL=C
export TZ=UTC
export RUST_LOG=error
export PARITY_RANDOM_SEED="${PARITY_RANDOM_SEED:-1337}"
export PYTHONHASHSEED="${PYTHONHASHSEED:-${PARITY_RANDOM_SEED}}"
export SOURCE_DATE_EPOCH="${SOURCE_DATE_EPOCH:-1700000000}"

run_java="${repo_root}/scripts/parity/run_java_gatk.sh"
run_rust="${repo_root}/scripts/parity/run_rust_gatk.sh"
diff_py="${repo_root}/scripts/parity/diff_outputs.py"
report_py="${repo_root}/scripts/parity/report.py"
smoke_profile="${PARITY_SMOKE_PROFILE:-smoke}"
if [[ "${smoke_profile}" != "smoke" && "${smoke_profile}" != "extended" ]]; then
  echo "Unsupported PARITY_SMOKE_PROFILE='${smoke_profile}' (expected: smoke|extended)" >&2
  exit 2
fi

check_labels=(
  "help-exit"
  "version-exit"
  "version-banner-content"
  "validate-help-mapped-content"
  "haplotypecaller-summary-content"
  "printreads-summary-content"
  "haplotypecaller-help-exit"
  "printreads-help-exit"
  "hc-invalid-interval-exit"
  "validate-sam-file-success"
  "validate-bam-file-success"
  "validate-vcf-file-success"
  "validate-fasta-file-success"
  "countbases-interval-chr1-1-16"
)
check_java_args=(
  "--help"
  "--version"
  "--version"
  "ValidateSamFile --help"
  "HaplotypeCaller --help"
  "PrintReads --help"
  "HaplotypeCaller --help"
  "PrintReads --help"
  "HaplotypeCaller -R parity/fixtures/reference.fa -I missing.bam -O out.vcf -L chr999:1-10"
  "ValidateSamFile -I parity/fixtures/sample.sam -MODE SUMMARY -IGNORE_WARNINGS true -REFERENCE_SEQUENCE parity/fixtures/reference.fa"
  "ValidateSamFile -I parity/fixtures/sample.bam -MODE SUMMARY -IGNORE_WARNINGS true -REFERENCE_SEQUENCE parity/fixtures/reference.fa"
  "ValidateVariants -V parity/fixtures/sample.vcf -R parity/fixtures/reference.fa"
  "CountBasesInReference -R parity/fixtures/reference.fa"
  "CountBasesInReference -R parity/fixtures/reference.fa -L chr1:1-16"
)
check_rust_args=(
  "--help"
  "--version"
  "--version"
  "Validate --help"
  "HaplotypeCaller --help"
  "PrintReads --help"
  "HaplotypeCaller --help"
  "PrintReads --help"
  "HaplotypeCaller -R parity/fixtures/reference.fa -I missing.bam -O out.vcf -L chr999:1-10"
  "Validate parity/fixtures/sample.sam -t SAM -R parity/fixtures/reference.fa"
  "Validate parity/fixtures/sample.bam -t BAM -R parity/fixtures/reference.fa"
  "Validate parity/fixtures/sample.vcf -t VCF -R parity/fixtures/reference.fa"
  "Validate parity/fixtures/reference.fa -t FASTA"
  "CountBasesInReference -R parity/fixtures/reference.fa -L chr1:1-16"
)
check_modes=(
  "exit-only"
  "exit-only"
  "normalized"
  "normalized"
  "normalized"
  "normalized"
  "exit-only"
  "exit-only"
  "exit-only"
  "normalized"
  "normalized"
  "normalized"
  "normalized"
  "normalized"
)
check_extract_regex=(
  ""
  ""
  "(?i)(The Genome Analysis Toolkit \\(GATK\\) v[0-9.]+|gatk-rs [0-9.]+ \\(independent community project)"
  "(?i)validat"
  "Call germline SNPs and indels via local re-assembly of haplotypes"
  "PrintReads"
  ""
  ""
  ""
  "(?i)(No errors found|SAM validation passed)"
  "(?i)(No errors found|BAM validation passed)"
  "(?i)(VCF validation passed|Processed [0-9]+ total variants)"
  "(?i)(FASTA validation passed|Processed [0-9]+ total bases)"
  "(?m)^[ACGTN]\\s*:\\s*[0-9]+"
)
check_presence_only=(
  "0"
  "0"
  "1"
  "1"
  "0"
  "1"
  "0"
  "0"
  "0"
  "1"
  "1"
  "1"
  "1"
  "0"
)
check_require_same_exit=(
  "1"
  "1"
  "1"
  "0"
  "1"
  "1"
  "1"
  "1"
  "1"
  "1"
  "1"
  "1"
  "1"
  "1"
)

if [[ "${smoke_profile}" == "extended" ]]; then
  check_labels+=(
    "countbases-interval-list-full"
    "countbases-whole-reference"
  )
  check_java_args+=(
    "CountBasesInReference -R parity/fixtures/reference.fa -L parity/fixtures/regions.interval_list"
    "CountBasesInReference -R parity/fixtures/reference.fa"
  )
  check_rust_args+=(
    "CountBasesInReference -R parity/fixtures/reference.fa -L parity/fixtures/regions.interval_list"
    "CountBasesInReference -R parity/fixtures/reference.fa"
  )
  check_modes+=(
    "normalized"
    "normalized"
  )
  check_extract_regex+=(
    "(?m)^[ACGTN]\\s*:\\s*[0-9]+"
    "(?m)^[ACGTN]\\s*:\\s*[0-9]+"
  )
  check_presence_only+=(
    "0"
    "0"
  )
  check_require_same_exit+=(
    "1"
    "1"
  )
fi

checks_json=()
passed=0
failed=0
skipped=0

for idx in "${!check_labels[@]}"; do
  label="${check_labels[$idx]}"
  java_args_string="${check_java_args[$idx]}"
  rust_args_string="${check_rust_args[$idx]}"
  mode="${check_modes[$idx]}"
  extract_regex="${check_extract_regex[$idx]}"
  presence_only="${check_presence_only[$idx]}"
  require_same_exit="${check_require_same_exit[$idx]}"
  IFS=' ' read -r -a java_args <<< "${java_args_string}"
  IFS=' ' read -r -a rust_args <<< "${rust_args_string}"

  java_out="${report_dir}/${label}.java.out"
  rust_out="${report_dir}/${label}.rust.out"
  check_json="${report_dir}/${label}.json"

  set +e
  "${run_java}" "${java_out}" "${java_args[@]}"
  java_exit=$?
  "${run_rust}" "${rust_out}" "${rust_args[@]}"
  rust_exit=$?
  set -e

  # Allow CI to skip differential checks when Java GATK is unavailable.
  if [[ "${java_exit}" -eq 127 && "${PARITY_ALLOW_MISSING_JAVA:-0}" == "1" ]]; then
    cat > "${check_json}" <<EOF
{
  "label": "${label}",
  "mode": "normalized",
  "equal": null,
  "skipped": true,
  "reason": "java_gatk_missing",
  "java_output": "${java_out}",
  "rust_output": "${rust_out}"
}
EOF
    skipped=$((skipped + 1))
    checks_json+=("${check_json}")
    continue
  fi

  # Exit code mismatch is considered a fail unless disabled for mapped checks.
  if [[ "${require_same_exit}" == "1" && "${java_exit}" -ne "${rust_exit}" ]]; then
    cat > "${check_json}" <<EOF
{
  "label": "${label}",
  "mode": "normalized",
  "equal": false,
  "reason": "exit_code_mismatch",
  "java_exit": ${java_exit},
  "rust_exit": ${rust_exit},
  "java_output": "${java_out}",
  "rust_output": "${rust_out}"
}
EOF
    failed=$((failed + 1))
    checks_json+=("${check_json}")
    echo "[parity-smoke] FAIL ${label} (exit_code_mismatch java=${java_exit} rust=${rust_exit})" >&2
    continue
  fi

  if [[ "${mode}" == "exit-only" ]]; then
    if [[ "${java_exit}" -eq "${rust_exit}" ]]; then
      cmp_exit=0
      cat > "${check_json}" <<EOF
{
  "label": "${label}",
  "mode": "${mode}",
  "equal": true,
  "java_exit": ${java_exit},
  "rust_exit": ${rust_exit},
  "java_output": "${java_out}",
  "rust_output": "${rust_out}"
}
EOF
    else
      cmp_exit=1
      cat > "${check_json}" <<EOF
{
  "label": "${label}",
  "mode": "${mode}",
  "equal": false,
  "reason": "exit_code_mismatch",
  "java_exit": ${java_exit},
  "rust_exit": ${rust_exit},
  "java_output": "${java_out}",
  "rust_output": "${rust_out}"
}
EOF
    fi
  else
    set +e
    diff_cmd=(
      python3 "${diff_py}"
      --java "${java_out}"
      --rust "${rust_out}"
      --label "${label}"
      --mode normalized
      --json-out "${check_json}"
    )
    if [[ -n "${extract_regex}" ]]; then
      diff_cmd+=("--extract-regex=${extract_regex}")
    fi
    if [[ "${presence_only}" == "1" ]]; then
      diff_cmd+=(--presence-only)
    fi
    "${diff_cmd[@]}"
    cmp_exit=$?
    set -e
  fi

  if [[ "${cmp_exit}" -eq 0 ]]; then
    passed=$((passed + 1))
    echo "[parity-smoke] ok ${label} (java_exit=${java_exit} rust_exit=${rust_exit})"
  else
    failed=$((failed + 1))
    echo "[parity-smoke] FAIL ${label} (java_exit=${java_exit} rust_exit=${rust_exit} cmp=${cmp_exit})" >&2
  fi

  checks_json+=("${check_json}")
done

# PrintReads: same input, compare output SAM (normalized; ignores @PG/@CO).
# Requires `RG` on reads so Java GATK does not filter them via WellformedReadFilter.
pr_label="printreads-sam-file-parity"
pr_java_stdout="${report_dir}/${pr_label}.java.stdout"
pr_rust_stdout="${report_dir}/${pr_label}.rust.stdout"
pr_java_sam="${report_dir}/${pr_label}.java.sam"
pr_rust_sam="${report_dir}/${pr_label}.rust.sam"
pr_check_json="${report_dir}/${pr_label}.json"
compare_sam_py="${repo_root}/scripts/parity/compare_sam_parity.py"

set +e
"${run_java}" "${pr_java_stdout}" PrintReads -I parity/fixtures/sample.sam -O "${pr_java_sam}"
pr_java_exit=$?
"${run_rust}" "${pr_rust_stdout}" PrintReads -I parity/fixtures/sample.sam -O "${pr_rust_sam}"
pr_rust_exit=$?
set -e

if [[ "${pr_java_exit}" -eq 127 && "${PARITY_ALLOW_MISSING_JAVA:-0}" == "1" ]]; then
  cat > "${pr_check_json}" <<EOF
{
  "label": "${pr_label}",
  "mode": "sam-file-parity",
  "equal": null,
  "skipped": true,
  "reason": "java_gatk_missing",
  "java_output": "${pr_java_stdout}",
  "rust_output": "${pr_rust_stdout}"
}
EOF
  skipped=$((skipped + 1))
  checks_json+=("${pr_check_json}")
elif [[ "${pr_java_exit}" -ne 0 || "${pr_rust_exit}" -ne 0 ]]; then
  cat > "${pr_check_json}" <<EOF
{
  "label": "${pr_label}",
  "mode": "sam-file-parity",
  "equal": false,
  "reason": "tool_exit_nonzero",
  "java_exit": ${pr_java_exit},
  "rust_exit": ${pr_rust_exit},
  "java_output": "${pr_java_stdout}",
  "rust_output": "${pr_rust_stdout}"
}
EOF
  failed=$((failed + 1))
  checks_json+=("${pr_check_json}")
  echo "[parity-smoke] FAIL ${pr_label} (tool_exit_nonzero java=${pr_java_exit} rust=${pr_rust_exit})" >&2
else
  set +e
  python3 "${compare_sam_py}" \
    --java-sam "${pr_java_sam}" \
    --rust-sam "${pr_rust_sam}" \
    --label "${pr_label}" \
    --json-out "${pr_check_json}"
  pr_cmp_exit=$?
  set -e
  if [[ "${pr_cmp_exit}" -eq 0 ]]; then
    passed=$((passed + 1))
    echo "[parity-smoke] ok ${pr_label}"
  else
    failed=$((failed + 1))
    echo "[parity-smoke] FAIL ${pr_label} (sam mismatch)" >&2
  fi
  checks_json+=("${pr_check_json}")
fi

# Strict VCF comparator self-check (Phase 0 / Step 7).
vcf_strict_label="vcf-strict-selfcheck"
vcf_strict_json="${report_dir}/${vcf_strict_label}.json"
strict_py="${repo_root}/scripts/parity/compare_vcf_strict.py"
set +e
python3 "${strict_py}" \
  --java "${repo_root}/parity/expected/sample.strict.vcf" \
  --rust "${repo_root}/parity/fixtures/sample.vcf" \
  --label "${vcf_strict_label}"
vs_exit=$?
set -e
if [[ "${vs_exit}" -eq 0 ]]; then
  passed=$((passed + 1))
  printf '%s\n' "{\"label\":\"${vcf_strict_label}\",\"mode\":\"vcf-strict\",\"equal\":true}" >"${vcf_strict_json}"
else
  failed=$((failed + 1))
  printf '%s\n' "{\"label\":\"${vcf_strict_label}\",\"mode\":\"vcf-strict\",\"equal\":false}" >"${vcf_strict_json}"
fi
checks_json+=("${vcf_strict_json}")

# Normalized VCF comparator self-check (Phase 0 / Step 8).
vcf_norm_label="vcf-normalized-selfcheck"
vcf_norm_json="${report_dir}/${vcf_norm_label}.json"
norm_py="${repo_root}/scripts/parity/compare_vcf_normalized.py"
set +e
python3 "${norm_py}" \
  --java "${repo_root}/parity/fixtures/sample_normalized_a.vcf" \
  --rust "${repo_root}/parity/fixtures/sample_normalized_b.vcf" \
  --label "${vcf_norm_label}"
vn_exit=$?
set -e
if [[ "${vn_exit}" -eq 0 ]]; then
  passed=$((passed + 1))
  printf '%s\n' "{\"label\":\"${vcf_norm_label}\",\"mode\":\"vcf-normalized\",\"equal\":true}" >"${vcf_norm_json}"
else
  failed=$((failed + 1))
  printf '%s\n' "{\"label\":\"${vcf_norm_label}\",\"mode\":\"vcf-normalized\",\"equal\":false}" >"${vcf_norm_json}"
fi
checks_json+=("${vcf_norm_json}")

# BAM alignment parity (Phase 0 / Step 9): requires `samtools` when PARITY_REQUIRE_SAMTOOLS=1 (CI).
bam_cmp_label="bam-alignment-selfcheck"
bam_cmp_json="${report_dir}/${bam_cmp_label}.json"
bam_cmp_py="${repo_root}/scripts/parity/compare_bam_alignment_parity.py"
set +e
python3 "${bam_cmp_py}" \
  --java-bam "${repo_root}/parity/fixtures/sample.bam" \
  --rust-bam "${repo_root}/parity/fixtures/sample.bam" \
  --label "${bam_cmp_label}" \
  --json-out "${bam_cmp_json}"
bam_exit=$?
set -e
if [[ "${bam_exit}" -eq 0 ]]; then
  passed=$((passed + 1))
  checks_json+=("${bam_cmp_json}")
elif [[ "${bam_exit}" -eq 1 ]]; then
  failed=$((failed + 1))
  checks_json+=("${bam_cmp_json}")
else
  failed=$((failed + 1))
  checks_json+=("${bam_cmp_json}")
fi

bam_cross_label="bam-alignment-crossformat-selfcheck"
bam_cross_json="${report_dir}/${bam_cross_label}.json"
set +e
python3 "${bam_cmp_py}" \
  --java-bam "${repo_root}/parity/fixtures/sample.bam" \
  --rust-bam "${repo_root}/parity/fixtures/sample.sam" \
  --label "${bam_cross_label}" \
  --json-out "${bam_cross_json}"
bam_cross_exit=$?
set -e
if [[ "${bam_cross_exit}" -eq 0 ]]; then
  passed=$((passed + 1))
else
  failed=$((failed + 1))
fi
checks_json+=("${bam_cross_json}")

# Step 24 locked expected snippets for reference+interval layer.
locked_countbases_check() {
  local label="$1"
  local expected_file="$2"
  local out_file="${report_dir}/${label}.rust.out"
  local check_json="${report_dir}/${label}.locked-expected.json"
  if [[ "${PARITY_ALLOW_MISSING_JAVA:-0}" == "1" ]]; then
    cat > "${check_json}" <<EOF
{
  "label": "${label}-locked-expected",
  "mode": "normalized",
  "equal": null,
  "skipped": true,
  "reason": "locked_expected_skipped_when_java_missing"
}
EOF
    skipped=$((skipped + 1))
    checks_json+=("${check_json}")
    return
  fi
  if [[ ! -f "${out_file}" ]]; then
    echo "[parity-smoke] FAIL ${label}-locked-expected (missing ${out_file})" >&2
    cat > "${check_json}" <<EOF
{
  "label": "${label}-locked-expected",
  "mode": "normalized",
  "equal": false,
  "reason": "missing_rust_output",
  "rust_output": "${out_file}"
}
EOF
    failed=$((failed + 1))
    checks_json+=("${check_json}")
    return
  fi
  set +e
  # Same extract as the live CountBases differential check — ignore tool banners.
  python3 "${diff_py}" \
    --java "${expected_file}" \
    --rust "${out_file}" \
    --label "${label}-locked-expected" \
    --mode normalized \
    --extract-regex='(?m)^[ACGTN]\s*:\s*[0-9]+' \
    --json-out "${check_json}"
  local cmp_exit=$?
  set -e
  if [[ "${cmp_exit}" -eq 0 ]]; then
    passed=$((passed + 1))
    echo "[parity-smoke] ok ${label}-locked-expected"
  else
    failed=$((failed + 1))
    echo "[parity-smoke] FAIL ${label}-locked-expected" >&2
  fi
  checks_json+=("${check_json}")
}

locked_countbases_check \
  "countbases-interval-chr1-1-16" \
  "${repo_root}/parity/expected/countbases-chr1-1-16.lines.txt"

if [[ "${smoke_profile}" == "extended" ]]; then
  locked_countbases_check \
    "countbases-interval-list-full" \
    "${repo_root}/parity/expected/countbases-interval-list-full.lines.txt"
  locked_countbases_check \
    "countbases-whole-reference" \
    "${repo_root}/parity/expected/countbases-whole-reference.lines.txt"
fi

summary_json="${report_dir}/parity-smoke.json"
{
  echo "{"
  echo "  \"profile\": \"${smoke_profile}\","
  echo "  \"determinism\": {\"rayon_threads\": \"${RAYON_NUM_THREADS}\", \"seed\": \"${PARITY_RANDOM_SEED}\", \"locale\": \"${LC_ALL}\", \"tz\": \"${TZ}\"},"
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

python3 "${report_py}" \
  --input-json "${summary_json}" \
  --output-md "${report_dir}/parity-smoke.md"

if [[ "${failed}" -gt 0 ]]; then
  echo "Parity smoke failed: ${failed} checks failed"
  exit 1
fi

echo "Parity smoke completed: ${passed} passed, ${skipped} skipped"
