#!/usr/bin/env bash
# Phase D.2.7 gate: contaminationFractionToFilter on isActive pileups (raw-activity-contam).
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

if [[ "${PARITY_SKIP_HC_FULL_D2C_CONTAM:-0}" == "1" ]]; then
  echo "[hc-full-parity-d2c-contam] skipped (PARITY_SKIP_HC_FULL_D2C_CONTAM=1)"
  exit 0
fi

cases_tsv="${repo_root}/parity/fixtures/hc-full-parity/d2c-contam/cases.tsv"
tmp_dir="${repo_root}/parity/reports/hc-full-parity-d2c-contam-tmp"
mkdir -p "${tmp_dir}"
cargo_run=(cargo run -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)

while IFS=$'\t' read -r case_id ref bam interval contam expected; do
  [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
  actual="${tmp_dir}/${case_id}.actual.tsv"
  "${cargo_run[@]}" raw-activity-contam \
    "${repo_root}/${ref}" "${repo_root}/${bam}" "${interval}" "${contam}" >"${actual}" 2>/dev/null
  if ! cmp -s "${repo_root}/${expected}" "${actual}"; then
    echo "[hc-full-parity-d2c-contam] mismatch case=${case_id}" >&2
    diff -u "${repo_root}/${expected}" "${actual}" >&2 || true
    exit 1
  fi
  echo "[hc-full-parity-d2c-contam] ok ${case_id}"
done <"${cases_tsv}"

echo "[hc-full-parity-d2c-contam] all cases match goldens."
