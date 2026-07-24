#!/usr/bin/env bash
# G-D05 — AlleleSubsettingUtils PL/AD/SAC/GQ vs frozen goldens.
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

cases_tsv="${repo_root}/parity/fixtures/hc-full-parity/g-subset-vc/cases.tsv"
tmp_dir="${repo_root}/parity/reports/hc-full-parity-g-subset-vc-tmp"
mkdir -p "${tmp_dir}"

cargo_run=(cargo run -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)
[[ "${PARITY_RUST_PROFILE:-dev}" == "release" ]] && cargo_run=(cargo run --release -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)

while IFS=$'\t' read -r case_id _pl _ad _sac expected; do
  [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
  actual="${tmp_dir}/${case_id}.actual.tsv"
  "${cargo_run[@]}" subset-alleles-vc \
    "${repo_root}/parity/fixtures/hc-full-parity/g-subset-vc/het_ac_sac_vc.tsv" >"${actual}" 2>"${tmp_dir}/${case_id}.stderr"
  cmp -s "${repo_root}/${expected}" "${actual}" || {
    echo "[hc-full-parity-g-subset-vc] mismatch ${case_id}" >&2
    diff -u "${repo_root}/${expected}" "${actual}" >&2 || true
    exit 1
  }
  echo "[hc-full-parity-g-subset-vc] ok ${case_id}"
done <"${cases_tsv}"

echo "[hc-full-parity-g-subset-vc] all cases match goldens."
