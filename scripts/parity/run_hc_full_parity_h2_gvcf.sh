#!/usr/bin/env bash
# Phase H.2.1 — GVCF header + writer block dump vs frozen goldens.
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

tmp_dir="${repo_root}/parity/reports/hc-full-parity-h2-tmp"
mkdir -p "${tmp_dir}"

cargo_run=(cargo run -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)
[[ "${PARITY_RUST_PROFILE:-dev}" == "release" ]] && cargo_run=(cargo run --release -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)

cases_tsv="${repo_root}/parity/fixtures/hc-full-parity/h2/cases.tsv"
while IFS=$'\t' read -r case_id contig length expected; do
  [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
  actual="${tmp_dir}/${case_id}.actual.tsv"
  "${cargo_run[@]}" gvcf-header "${contig}" "${length}" >"${actual}"
  cmp -s "${repo_root}/${expected}" "${actual}" || {
    echo "[hc-full-parity-h2] header mismatch ${case_id}" >&2
    exit 1
  }
  echo "[hc-full-parity-h2] ok header ${case_id}"
done <"${cases_tsv}"

cases_tsv="${repo_root}/parity/fixtures/hc-full-parity/h2-blocks/cases.tsv"
while IFS=$'\t' read -r case_id fixture expected; do
  [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
  actual="${tmp_dir}/${case_id}_blocks.actual.tsv"
  "${cargo_run[@]}" gvcf-writer-blocks "${repo_root}/${fixture}" >"${actual}"
  cmp -s "${repo_root}/${expected}" "${actual}" || {
    echo "[hc-full-parity-h2] blocks mismatch ${case_id}" >&2
    diff -u "${repo_root}/${expected}" "${actual}" >&2 || true
    exit 1
  }
  echo "[hc-full-parity-h2] ok blocks ${case_id}"
done <"${cases_tsv}"

echo "[hc-full-parity-h2] all cases match goldens."
