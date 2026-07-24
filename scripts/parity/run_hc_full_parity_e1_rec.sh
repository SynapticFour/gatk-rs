#!/usr/bin/env bash
# Phase E.1.1: PileupReadErrorCorrector vs frozen goldens.
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

if [[ "${PARITY_SKIP_HC_FULL_E1_REC:-0}" == "1" ]]; then
  echo "[hc-full-parity-e1-rec] skipped (PARITY_SKIP_HC_FULL_E1_REC=1)"
  exit 0
fi

cases_tsv="${repo_root}/parity/fixtures/hc-full-parity/e1-rec/cases.tsv"
tmp_dir="${repo_root}/parity/reports/hc-full-parity-e1-rec-tmp"
mkdir -p "${tmp_dir}"

cargo_run=(cargo run -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)
if [[ "${PARITY_RUST_PROFILE:-dev}" == "release" ]]; then
  cargo_run=(cargo run --release -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)
fi

while IFS=$'\t' read -r case_id reads log_odds expected; do
  [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
  actual="${tmp_dir}/${case_id}.actual.tsv"
  "${cargo_run[@]}" read-error-correction \
    "${repo_root}/${reads}" \
    "${log_odds}" \
    >"${actual}" 2>"${tmp_dir}/${case_id}.stderr"
  if ! cmp -s "${repo_root}/${expected}" "${actual}"; then
    echo "[hc-full-parity-e1-rec] mismatch case=${case_id}" >&2
    diff -u "${repo_root}/${expected}" "${actual}" >&2 || true
    exit 1
  fi
  echo "[hc-full-parity-e1-rec] ok ${case_id}"
done <"${cases_tsv}"

echo "[hc-full-parity-e1-rec] all cases match goldens."
