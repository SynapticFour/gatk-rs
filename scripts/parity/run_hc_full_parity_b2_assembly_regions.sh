#!/usr/bin/env bash
# Phase B.2 gate: Rust assembly-region iterator TSV matches frozen goldens (parity/fixtures/hc-full-parity/b2/README.md).
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

if [[ "${PARITY_SKIP_HC_FULL_B2:-0}" == "1" ]]; then
  echo "[hc-full-parity-b2] skipped (PARITY_SKIP_HC_FULL_B2=1)"
  exit 0
fi

cases_tsv="${repo_root}/parity/fixtures/hc-full-parity/b2/cases.tsv"
report_dir="${repo_root}/parity/reports"
tmp_dir="${report_dir}/hc-full-parity-b2-tmp"
mkdir -p "${tmp_dir}"

profile="${PARITY_RUST_PROFILE:-dev}"
cargo_run=(cargo run -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)
if [[ "${profile}" == "release" ]]; then
  cargo_run=(cargo run --release -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)
fi

while IFS=$'\t' read -r case_id ref bam interval padding expected; do
  [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
  actual="${tmp_dir}/${case_id}.actual.tsv"
  set +e
  if [[ -z "${padding}" || "${padding}" == "-" ]]; then
    "${cargo_run[@]}" assembly-regions "${repo_root}/${ref}" "${repo_root}/${bam}" "${interval}" >"${actual}" 2>"${tmp_dir}/${case_id}.stderr"
  else
    "${cargo_run[@]}" assembly-regions "${repo_root}/${ref}" "${repo_root}/${bam}" "${interval}" "${padding}" >"${actual}" 2>"${tmp_dir}/${case_id}.stderr"
  fi
  dump_ec=$?
  set -e
  if [[ "${dump_ec}" -ne 0 ]]; then
    echo "[hc-full-parity-b2] dump failed (case=${case_id}, exit=${dump_ec}). stderr:" >&2
    cat "${tmp_dir}/${case_id}.stderr" >&2
    exit 1
  fi
  if ! cmp -s "${repo_root}/${expected}" "${actual}"; then
    echo "[hc-full-parity-b2] mismatch case=${case_id}" >&2
    echo "  expected: ${repo_root}/${expected}" >&2
    echo "  actual:   ${actual}" >&2
    diff -u "${repo_root}/${expected}" "${actual}" >&2 || true
    exit 1
  fi
  echo "[hc-full-parity-b2] ok ${case_id}"
done <"${cases_tsv}"

echo "[hc-full-parity-b2] all cases match goldens."
