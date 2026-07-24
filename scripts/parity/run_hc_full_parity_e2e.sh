#!/usr/bin/env bash
# Phase E2E.1: assembly region → haplotypes (active/inactive region from BAM).
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

if [[ "${PARITY_SKIP_HC_FULL_E2E:-0}" == "1" ]]; then
  echo "[hc-full-parity-e2e] skipped (PARITY_SKIP_HC_FULL_E2E=1)"
  exit 0
fi

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

cases_tsv="${repo_root}/parity/fixtures/hc-full-parity/e2e/cases.tsv"
tmp_dir="${repo_root}/parity/reports/hc-full-parity-e2e-tmp"
mkdir -p "${tmp_dir}"

cargo_run=(cargo run -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)
if [[ "${PARITY_RUST_PROFILE:-dev}" == "release" ]]; then
  cargo_run=(cargo run --release -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)
fi

while IFS=$'\t' read -r case_id ref bam interval padding target _l2_strict expected; do
  [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
  actual="${tmp_dir}/${case_id}.actual.tsv"
  local_fa="${repo_root}/${ref}"
  if [[ -f "${local_fa}" && ! -f "${local_fa}.fai" ]] && command -v samtools >/dev/null 2>&1; then
    samtools faidx "${local_fa}" 2>/dev/null || true
  fi
  bam_resolved="$(resolve_alignment_path "${repo_root}/${bam}")"
  args=(
    assembly-region-haplotypes
    "${local_fa}"
    "${bam_resolved}"
    "${interval}"
  )
  [[ -n "${padding}" && "${padding}" != "-" ]] && args+=("${padding}")
  [[ -n "${target}" && "${target}" != "-" && "${target}" != "active" ]] && args+=("${target}")
  "${cargo_run[@]}" "${args[@]}" >"${actual}" 2>"${tmp_dir}/${case_id}.stderr"
  if ! cmp -s "${repo_root}/${expected}" "${actual}"; then
    echo "[hc-full-parity-e2e] mismatch case=${case_id}" >&2
    diff -u "${repo_root}/${expected}" "${actual}" >&2 || true
    exit 1
  fi
  echo "[hc-full-parity-e2e] ok ${case_id}"
done <"${cases_tsv}"

echo "[hc-full-parity-e2e] all cases match goldens."
