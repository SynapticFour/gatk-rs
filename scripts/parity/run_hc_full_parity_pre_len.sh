#!/usr/bin/env bash
# Phase PRE.2 gate: assembly-path read length filter (unclippedReadLength >= 10).
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

if [[ "${PARITY_SKIP_HC_FULL_PRE_LEN:-0}" == "1" ]]; then
  echo "[hc-full-parity-pre-len] skipped (PARITY_SKIP_HC_FULL_PRE_LEN=1)"
  exit 0
fi

cases_tsv="${repo_root}/parity/fixtures/hc-full-parity/pre-len/cases.tsv"
tmp_dir="${repo_root}/parity/reports/hc-full-parity-pre-len-tmp"
mkdir -p "${tmp_dir}"
cargo_run=(cargo run -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)

while IFS=$'\t' read -r case_id bam expected; do
  [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
  actual="${tmp_dir}/${case_id}.actual.tsv"
  "${cargo_run[@]}" read-pre-len "${repo_root}/${bam}" >"${actual}" 2>/dev/null
  if ! cmp -s "${repo_root}/${expected}" "${actual}"; then
    echo "[hc-full-parity-pre-len] mismatch case=${case_id}" >&2
    diff -u "${repo_root}/${expected}" "${actual}" >&2 || true
    exit 1
  fi
  echo "[hc-full-parity-pre-len] ok ${case_id}"
done <"${cases_tsv}"

echo "[hc-full-parity-pre-len] all cases match goldens."
