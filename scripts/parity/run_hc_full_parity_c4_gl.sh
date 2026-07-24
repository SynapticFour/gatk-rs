#!/usr/bin/env bash
# Phase C.4 — per-locus genotype likelihoods + activity (GAP-C-01 L2).
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

if [[ "${PARITY_SKIP_HC_FULL_C4_GL:-0}" == "1" ]]; then
  echo "[hc-full-parity-c4-gl] skipped (PARITY_SKIP_HC_FULL_C4_GL=1)"
  exit 0
fi

cases_tsv="${repo_root}/parity/fixtures/hc-full-parity/c4-gl/cases.tsv"
report_dir="${repo_root}/parity/reports"
tmp_dir="${report_dir}/hc-full-parity-c4-gl-tmp"
mkdir -p "${tmp_dir}"

profile="${PARITY_RUST_PROFILE:-dev}"
cargo_run=(cargo run -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)
if [[ "${profile}" == "release" ]]; then
  cargo_run=(cargo run --release -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)
fi

while IFS=$'\t' read -r case_id ref bam interval _expected; do
  [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
  actual="${tmp_dir}/${case_id}.actual.tsv"
  set +e
  "${cargo_run[@]}" genotype-likelihoods \
    "${repo_root}/${ref}" "${repo_root}/${bam}" "${interval}" "-" >"${actual}" 2>"${tmp_dir}/${case_id}.stderr"
  dump_ec=$?
  set -e
  if [[ "${dump_ec}" -ne 0 ]]; then
    echo "[hc-full-parity-c4-gl] dump failed (case=${case_id})" >&2
    cat "${tmp_dir}/${case_id}.stderr" >&2
    exit 1
  fi
  if ! cmp -s "${repo_root}/${_expected}" "${actual}"; then
    echo "[hc-full-parity-c4-gl] mismatch case=${case_id}" >&2
    diff -u "${repo_root}/${_expected}" "${actual}" >&2 || true
    exit 1
  fi
  echo "[hc-full-parity-c4-gl] ok ${case_id}"
done <"${cases_tsv}"

echo "[hc-full-parity-c4-gl] all cases match goldens."
