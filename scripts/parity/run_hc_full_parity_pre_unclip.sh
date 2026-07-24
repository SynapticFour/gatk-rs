#!/usr/bin/env bash
# Phase PRE.1 gate: HC soft-clip policy (revert vs hard-clip) before assembly.
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

if [[ "${PARITY_SKIP_HC_FULL_PRE_UNCLIP:-0}" == "1" ]]; then
  echo "[hc-full-parity-pre-unclip] skipped (PARITY_SKIP_HC_FULL_PRE_UNCLIP=1)"
  exit 0
fi

cases_tsv="${repo_root}/parity/fixtures/hc-full-parity/pre/cases.tsv"
tmp_dir="${repo_root}/parity/reports/hc-full-parity-pre-tmp"
mkdir -p "${tmp_dir}"
cargo_run=(cargo run -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)

while IFS=$'\t' read -r case_id bam expected; do
  [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
  actual="${tmp_dir}/${case_id}.actual.tsv"
  "${cargo_run[@]}" read-pre-softclip "${repo_root}/${bam}" 0 0 >"${actual}" 2>/dev/null
  if ! cmp -s "${repo_root}/${expected}" "${actual}"; then
    echo "[hc-full-parity-pre-unclip] mismatch case=${case_id}" >&2
    diff -u "${repo_root}/${expected}" "${actual}" >&2 || true
    exit 1
  fi
  echo "[hc-full-parity-pre-unclip] ok ${case_id}"
done <"${cases_tsv}"

echo "[hc-full-parity-pre-unclip] all cases match goldens."
