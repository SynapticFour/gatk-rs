#!/usr/bin/env bash
# Phase D.2 gate: positional + allele-biased downsampling summaries.
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

if [[ "${PARITY_SKIP_HC_FULL_D2:-0}" == "1" ]]; then
  echo "[hc-full-parity-d2] skipped (PARITY_SKIP_HC_FULL_D2=1)"
  exit 0
fi

cases_tsv="${repo_root}/parity/fixtures/hc-full-parity/d2/cases.tsv"
tmp_dir="${repo_root}/parity/reports/hc-full-parity-d2-tmp"
mkdir -p "${tmp_dir}"
cargo_run=(cargo run -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)

while IFS=$'\t' read -r case_id ref bam contig pos1 cap expected mode; do
  [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
  actual="${tmp_dir}/${case_id}.actual.tsv"
  if [[ "${case_id}" == positional_* ]]; then
    if [[ -n "${mode:-}" ]]; then
      "${cargo_run[@]}" downsample-positional "${repo_root}/${bam}" "${cap}" "${mode}" >"${actual}"
    else
      "${cargo_run[@]}" downsample-positional "${repo_root}/${bam}" "${cap}" >"${actual}"
    fi
  else
    "${cargo_run[@]}" downsample-allele "${repo_root}/${ref}" "${repo_root}/${bam}" "${contig}" "${pos1}" "${cap}" >"${actual}"
  fi
  if ! cmp -s "${repo_root}/${expected}" "${actual}"; then
    echo "[hc-full-parity-d2] mismatch case=${case_id}" >&2
    diff -u "${repo_root}/${expected}" "${actual}" >&2 || true
    exit 1
  fi
  echo "[hc-full-parity-d2] ok ${case_id}"
done <"${cases_tsv}"

echo "[hc-full-parity-d2] all cases match goldens."
