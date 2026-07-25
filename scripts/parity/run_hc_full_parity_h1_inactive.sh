#!/usr/bin/env bash
# Phase H.1.2 — inactive region referenceModelForNoVariation vs frozen goldens.
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

cases_tsv="${repo_root}/parity/fixtures/hc-full-parity/h1-inactive/cases.tsv"
tmp_dir="${repo_root}/parity/reports/hc-full-parity-h1-inactive-tmp"
mkdir -p "${tmp_dir}"

cargo_run=(cargo run -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)
[[ "${PARITY_RUST_PROFILE:-dev}" == "release" ]] && cargo_run=(cargo run --release -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)

ran=0
skipped=0
while IFS=$'\t' read -r case_id ref bam interval padding expected; do
  [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
  if [[ ! -f "${repo_root}/${ref}" || ! -f "${repo_root}/${bam}" ]]; then
    echo "[hc-full-parity-h1-inactive] skip ${case_id}: missing ref/bam (realworld assets not staged)"
    skipped=$((skipped + 1))
    continue
  fi
  actual="${tmp_dir}/${case_id}.actual.tsv"
  "${cargo_run[@]}" inactive-reference-model \
    "${repo_root}/${ref}" "${repo_root}/${bam}" "${interval}" "${padding}" >"${actual}"
  cmp -s "${repo_root}/${expected}" "${actual}" || {
    echo "[hc-full-parity-h1-inactive] mismatch ${case_id}" >&2
    diff -u "${repo_root}/${expected}" "${actual}" >&2 || true
    exit 1
  }
  echo "[hc-full-parity-h1-inactive] ok ${case_id}"
  ran=$((ran + 1))
done <"${cases_tsv}"

if [[ "${ran}" -eq 0 ]]; then
  echo "[hc-full-parity-h1-inactive] ERROR: no cases ran (skipped=${skipped})" >&2
  exit 1
fi
echo "[hc-full-parity-h1-inactive] all ran cases match goldens (ran=${ran} skipped=${skipped})."
