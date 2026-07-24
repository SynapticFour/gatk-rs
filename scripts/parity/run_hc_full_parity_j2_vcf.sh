#!/usr/bin/env bash
# Phase J.2.3 — VCF identity fields (CHROM/POS/REF/ALT/QUAL/FILTER) parity dumps.
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

resolve_alignment_path() {
  local path="$1"
  if [[ "${path}" != *.sam ]]; then
    printf '%s\n' "${path}"
    return
  fi
  local cache_dir="${repo_root}/parity/build/sam-indexed-bam"
  mkdir -p "${cache_dir}"
  local out="${cache_dir}/$(basename "${path%.sam}").bam"
  if [[ ! -f "${out}" ]]; then
    samtools view -bS "${path}" | samtools sort -o "${out}"
    samtools index "${out}"
  fi
  printf '%s\n' "${out}"
}

cases_tsv="${repo_root}/parity/fixtures/hc-full-parity/j2/cases_vcf.tsv"
tmp_dir="${repo_root}/parity/reports/hc-full-parity-j2-tmp"
mkdir -p "${tmp_dir}"

cargo_run=(cargo run -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)
[[ "${PARITY_RUST_PROFILE:-dev}" == "release" ]] && cargo_run=(cargo run --release -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)

while IFS=$'\t' read -r case_id kind ref bam_or_contig interval_or_pos ref_allele alt_allele gl_csv ad_csv expected; do
  [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
  actual="${tmp_dir}/${case_id}.vcf.actual.tsv"
  if [[ "${kind}" == "call-region" ]]; then
    bam_resolved="$(resolve_alignment_path "${repo_root}/${bam_or_contig}")"
    "${cargo_run[@]}" call-region-vcf \
      "${repo_root}/${ref}" "${bam_resolved}" "${interval_or_pos}" >"${actual}"
  else
    "${cargo_run[@]}" variant-vcf-from-gl-ad \
      "${bam_or_contig}" "${interval_or_pos}" "${ref_allele}" "${alt_allele}" \
      "${gl_csv}" "${ad_csv}" >"${actual}"
  fi
  cmp -s "${repo_root}/${expected}" "${actual}" || {
    echo "[hc-full-parity-j2-vcf] mismatch ${case_id}" >&2
    diff -u "${repo_root}/${expected}" "${actual}" >&2 || true
    exit 1
  }
  echo "[hc-full-parity-j2-vcf] ok ${case_id}"
done <"${cases_tsv}"

echo "[hc-full-parity-j2-vcf] all cases match goldens."
