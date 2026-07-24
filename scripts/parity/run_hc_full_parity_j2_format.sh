#!/usr/bin/env bash
# Phase J.2.4 — agreed FORMAT subset (GT/GQ/DP/AD/PL) parity dumps.
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

cases_tsv="${repo_root}/parity/fixtures/hc-full-parity/j2/cases_format.tsv"
tmp_dir="${repo_root}/parity/reports/hc-full-parity-j2-tmp"
mkdir -p "${tmp_dir}"

cargo_run=(cargo run -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)
[[ "${PARITY_RUST_PROFILE:-dev}" == "release" ]] && cargo_run=(cargo run --release -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)

while IFS=$'\t' read -r case_id contig pos ref_allele alt_allele gl_csv ad_csv expected; do
  [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
  actual="${tmp_dir}/${case_id}.format.actual.tsv"
  "${cargo_run[@]}" variant-format-from-gl-ad \
    "${contig}" "${pos}" "${ref_allele}" "${alt_allele}" "${gl_csv}" "${ad_csv}" >"${actual}"
  cmp -s "${repo_root}/${expected}" "${actual}" || {
    echo "[hc-full-parity-j2-format] mismatch ${case_id}" >&2
    diff -u "${repo_root}/${expected}" "${actual}" >&2 || true
    exit 1
  }
  echo "[hc-full-parity-j2-format] ok ${case_id}"
done <"${cases_tsv}"

echo "[hc-full-parity-j2-format] all cases match goldens."
