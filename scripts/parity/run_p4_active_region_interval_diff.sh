#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
report_dir="${repo_root}/parity/reports"
mkdir -p "${report_dir}"
tmp_dir="${report_dir}/p4-assembly-region-diff-tmp"
mkdir -p "${tmp_dir}"

run_java="${repo_root}/scripts/parity/run_java_gatk.sh"
normalize_py="${repo_root}/scripts/parity/normalize_assembly_region_igv.py"
cases_tsv="${repo_root}/parity/fixtures/p4_assembly_region_cases.tsv"
ref_fa="${repo_root}/parity/fixtures/reference.fa"
in_bam="${repo_root}/parity/fixtures/sample.bam"

if [[ ! -f "${in_bam}.bai" ]]; then
  if command -v samtools >/dev/null 2>&1; then
    samtools index "${in_bam}"
  else
    echo "Missing ${in_bam}.bai and samtools on PATH; cannot index fixture BAM." >&2
    exit 2
  fi
fi

while IFS=$'\t' read -r case_id interval expected_rel; do
  [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
  expected="${repo_root}/${expected_rel}"
  java_stdout="${tmp_dir}/haplotypecaller_assembly_region_${case_id}.java.stdout"
  java_igv="${tmp_dir}/assembly_regions_${case_id}.java.igv"
  java_norm="${tmp_dir}/assembly_regions_${case_id}.java.normalized.txt"
  out_vcf="${tmp_dir}/hc_p4_scratch_${case_id}.vcf"

  set +e
  "${run_java}" "${java_stdout}" HaplotypeCaller \
    -R "${ref_fa}" \
    -I "${in_bam}" \
    -O "${out_vcf}" \
    -L "${interval}" \
    --assembly-region-out "${java_igv}"
  java_exit=$?
  set -e

  if [[ "${java_exit}" -eq 127 && "${PARITY_ALLOW_MISSING_JAVA:-0}" == "1" ]]; then
    echo "[p4-assembly-region-diff] skipped: Java GATK not available (exit 127, PARITY_ALLOW_MISSING_JAVA=1)"
    exit 0
  fi

  if [[ "${java_exit}" -ne 0 ]]; then
    echo "Java HaplotypeCaller failed (exit=${java_exit}, case=${case_id}). See ${java_stdout}" >&2
    exit 1
  fi

  if [[ ! -f "${java_igv}" ]]; then
    echo "Expected IGV output at ${java_igv} missing (case=${case_id})." >&2
    exit 1
  fi

  python3 "${normalize_py}" "${java_igv}" -o "${java_norm}"

  if ! cmp -s "${java_norm}" "${expected}"; then
    echo "Assembly region IGV normalized output differs from expected (case=${case_id})." >&2
    echo "  actual:   ${java_norm}" >&2
    echo "  expected: ${expected}" >&2
    diff -u "${expected}" "${java_norm}" >&2 || true
    exit 1
  fi
done < "${cases_tsv}"

echo "[p4-assembly-region-diff] Java HaplotypeCaller --assembly-region-out matches frozen expected corpus."
