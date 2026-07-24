#!/usr/bin/env bash
# Phase D.4 gate: HC shard pipeline — IUPAC pre-transform, read filter, post-transform (identity).
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

if [[ "${PARITY_SKIP_HC_FULL_D4:-0}" == "1" ]]; then
  echo "[hc-full-parity-d4] skipped (PARITY_SKIP_HC_FULL_D4=1)"
  exit 0
fi

cases_tsv="${repo_root}/parity/fixtures/hc-full-parity/d4/cases.tsv"
tmp_dir="${repo_root}/parity/reports/hc-full-parity-d4-tmp"
mkdir -p "${tmp_dir}"
cargo_run=(cargo run -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)

while IFS=$'\t' read -r case_id bam expected; do
  [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
  actual="${tmp_dir}/${case_id}.actual.tsv"
  "${cargo_run[@]}" read-shard-pipeline "${repo_root}/${bam}" >"${actual}" 2>/dev/null
  if ! cmp -s "${repo_root}/${expected}" "${actual}"; then
    echo "[hc-full-parity-d4] mismatch case=${case_id}" >&2
    diff -u "${repo_root}/${expected}" "${actual}" >&2 || true
    exit 1
  fi
  echo "[hc-full-parity-d4] ok ${case_id}"
done <"${cases_tsv}"

echo "[hc-full-parity-d4] all cases match goldens."
