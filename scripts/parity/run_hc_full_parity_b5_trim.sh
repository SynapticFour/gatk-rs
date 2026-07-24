#!/usr/bin/env bash
# Phase B.5.5 — AssemblyRegionTrimmer (post-assembly genotyping span).
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

if [[ "${PARITY_SKIP_HC_FULL_B5_TRIM:-0}" == "1" ]]; then
  echo "[hc-full-parity-b5-trim] skipped (PARITY_SKIP_HC_FULL_B5_TRIM=1)"
  exit 0
fi

cases_tsv="${repo_root}/parity/fixtures/hc-full-parity/b5-trim/cases.tsv"
report_dir="${repo_root}/parity/reports"
tmp_dir="${report_dir}/hc-full-parity-b5-trim-tmp"
mkdir -p "${tmp_dir}"

profile="${PARITY_RUST_PROFILE:-dev}"
cargo_run=(cargo run -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)
if [[ "${profile}" == "release" ]]; then
  cargo_run=(cargo run --release -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)
fi

while IFS=$'\t' read -r case_id ref contig start end ext_start ext_end variants legacy expected; do
  [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
  actual="${tmp_dir}/${case_id}.actual.tsv"
  args=(
    assembly-region-trim
    "${repo_root}/${ref}"
    "${contig}" "${start}" "${end}" "${ext_start}" "${ext_end}"
  )
  if [[ -n "${variants}" && "${variants}" != "-" ]]; then
    args+=("${repo_root}/${variants}")
  else
    args+=("-")
  fi
  if [[ "${legacy}" == "1" ]]; then
    args+=("legacy")
  fi
  set +e
  "${cargo_run[@]}" "${args[@]}" >"${actual}" 2>"${tmp_dir}/${case_id}.stderr"
  dump_ec=$?
  set -e
  if [[ "${dump_ec}" -ne 0 ]]; then
    echo "[hc-full-parity-b5-trim] dump failed (case=${case_id})" >&2
    cat "${tmp_dir}/${case_id}.stderr" >&2
    exit 1
  fi
  if ! cmp -s "${repo_root}/${expected}" "${actual}"; then
    echo "[hc-full-parity-b5-trim] mismatch case=${case_id}" >&2
    diff -u "${repo_root}/${expected}" "${actual}" >&2 || true
    exit 1
  fi
  echo "[hc-full-parity-b5-trim] ok ${case_id}"
done <"${cases_tsv}"

echo "[hc-full-parity-b5-trim] all cases match goldens."
