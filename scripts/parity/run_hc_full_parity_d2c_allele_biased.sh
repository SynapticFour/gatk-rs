#!/usr/bin/env bash
# Phase D.2.6 gate: GATK AlleleBiasedDownsamplingUtils (target counts + evidence removal).
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

if [[ "${PARITY_SKIP_HC_FULL_D2C:-0}" == "1" ]]; then
  echo "[hc-full-parity-d2c] skipped (PARITY_SKIP_HC_FULL_D2C=1)"
  exit 0
fi

cases_tsv="${repo_root}/parity/fixtures/hc-full-parity/d2c/cases.tsv"
tmp_dir="${repo_root}/parity/reports/hc-full-parity-d2c-tmp"
mkdir -p "${tmp_dir}"
cargo_run=(cargo run -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)

while IFS=$'\t' read -r case_id col2 col3 col4 col5 col6 col7; do
  [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
  actual="${tmp_dir}/${case_id}.actual.tsv"
  if [[ "${case_id}" == target_* ]]; then
    expected="${repo_root}/${col4}"
    "${cargo_run[@]}" allele-biased-target-counts "${col2}" "${col3}" >"${actual}" 2>/dev/null
  else
    expected="${repo_root}/${col7}"
    "${cargo_run[@]}" allele-biased-evidence \
      "${repo_root}/${col2}" "${repo_root}/${col3}" "${col4}" "${col5}" "${col6}" >"${actual}" 2>/dev/null
  fi
  if ! cmp -s "${expected}" "${actual}"; then
    echo "[hc-full-parity-d2c] mismatch case=${case_id}" >&2
    diff -u "${expected}" "${actual}" >&2 || true
    exit 1
  fi
  echo "[hc-full-parity-d2c] ok ${case_id}"
done <"${cases_tsv}"

echo "[hc-full-parity-d2c] all cases match goldens."
