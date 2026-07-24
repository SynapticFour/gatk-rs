#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
report_dir="${repo_root}/parity/reports"
mkdir -p "${report_dir}"
tmp_dir="${report_dir}/p3-region-records-diff-tmp"
mkdir -p "${tmp_dir}"

run_java="${repo_root}/scripts/parity/run_java_gatk.sh"
run_rust="${repo_root}/scripts/parity/run_rust_gatk.sh"

source_sam="${repo_root}/parity/fixtures/p3_region_reads.sam"
source_bam="${tmp_dir}/p3_region_reads.bam"
samtools view -bS "${source_sam}" > "${source_bam}"
samtools index "${source_bam}"

extract_java_qnames() {
  local path="$1"
  python3 - <<PY
import pathlib
out = []
for ln in pathlib.Path("${path}").read_text(encoding="utf-8", errors="replace").splitlines():
    if not ln or ln.startswith("@"):
        continue
    out.append(ln.split("\t", 1)[0])
print("\n".join(out))
PY
}

extract_rust_qnames() {
  local path="$1"
  python3 - <<PY
import pathlib
lines = [ln.strip() for ln in pathlib.Path("${path}").read_text(encoding="utf-8", errors="replace").splitlines() if ln.strip()]
print("\n".join(lines))
PY
}

passed=0
failed=0
skipped=0
checks_json=()

run_case() {
  local label="$1"
  local region="$2"
  local java_stdout="${report_dir}/${label}.java.stdout"
  local rust_stdout="${report_dir}/${label}.rust.stdout"
  local java_sam="${report_dir}/${label}.java.sam"
  local check_json="${report_dir}/${label}.json"

  set +e
  "${run_java}" "${java_stdout}" PrintReads -I "${source_bam}" -L "${region}" -O "${java_sam}"
  java_exit=$?
  "${run_rust}" "${rust_stdout}" ListReadsInRegion -I "${source_bam}" -L "${region}"
  rust_exit=$?
  set -e

  if [[ "${java_exit}" -eq 127 && "${PARITY_ALLOW_MISSING_JAVA:-0}" == "1" ]]; then
    cat > "${check_json}" <<EOF
{
  "label": "${label}",
  "mode": "phase3-region-records-diff",
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
  "mode": "phase3-region-records-diff",
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

  java_qnames="$(extract_java_qnames "${java_sam}")"
  rust_qnames="$(extract_rust_qnames "${rust_stdout}")"
  equal=false
  if [[ "${java_qnames}" == "${rust_qnames}" ]]; then
    equal=true
  fi

  cat > "${check_json}" <<EOF
{
  "label": "${label}",
  "mode": "phase3-region-records-diff",
  "equal": ${equal},
  "region": "${region}",
  "java_qnames": "$(echo "${java_qnames}" | tr '\n' ',' | sed 's/,$//')",
  "rust_qnames": "$(echo "${rust_qnames}" | tr '\n' ',' | sed 's/,$//')"
}
EOF

  if [[ "${equal}" == "true" ]]; then
    passed=$((passed + 1))
  else
    failed=$((failed + 1))
  fi
  checks_json+=("${check_json}")
}

run_case "p3-region-records-chr1-1-16" "chr1:1-16"
run_case "p3-region-records-chr1-17-32" "chr1:17-32"

summary_json="${report_dir}/p3-region-records-diff.json"
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
  echo "P3 region records differential checks failed: ${failed}"
  exit 1
fi

echo "P3 region records differential checks passed: ${passed}, skipped: ${skipped}"
