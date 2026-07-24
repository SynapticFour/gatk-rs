#!/usr/bin/env bash
# G-D05 — AlleleSubsettingUtils.subsetAlleles PL/AD vs frozen goldens.
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

cases_tsv="${repo_root}/parity/fixtures/hc-full-parity/g-subset-pl/cases.tsv"
tmp_dir="${repo_root}/parity/reports/hc-full-parity-g-subset-pl-tmp"
mkdir -p "${tmp_dir}"

cargo_run=(cargo run -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)
[[ "${PARITY_RUST_PROFILE:-dev}" == "release" ]] && cargo_run=(cargo run --release -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)

while IFS=$'\t' read -r case_id pl_csv ad_csv expected; do
  [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
  input="${tmp_dir}/${case_id}.input.tsv"
  printf '%s\t%s\t%s\n' "${case_id}" "${pl_csv}" "${ad_csv}" >"${input}"
  actual="${tmp_dir}/${case_id}.actual.tsv"
  "${cargo_run[@]}" subset-alleles-pl "${input}" >"${actual}" 2>"${tmp_dir}/${case_id}.stderr"
  cmp -s "${repo_root}/${expected}" "${actual}" || {
    echo "[hc-full-parity-g-subset-pl] mismatch ${case_id}" >&2
    diff -u "${repo_root}/${expected}" "${actual}" >&2 || true
    exit 1
  }
  echo "[hc-full-parity-g-subset-pl] ok ${case_id}"
done <"${cases_tsv}"

echo "[hc-full-parity-g-subset-pl] all cases match goldens."
