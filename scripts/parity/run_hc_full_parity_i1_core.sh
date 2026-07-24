#!/usr/bin/env bash
# Phase I.1.1 — core INFO annotations (AC/AN/AF/NS/DP) via VariantAnnotatorEngine parity v1.
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

cases_tsv="${repo_root}/parity/fixtures/hc-full-parity/i1/cases.tsv"
tmp_dir="${repo_root}/parity/reports/hc-full-parity-i1-tmp"
mkdir -p "${tmp_dir}"

cargo_run=(cargo run -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)
[[ "${PARITY_RUST_PROFILE:-dev}" == "release" ]] && cargo_run=(cargo run --release -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)

while IFS=$'\t' read -r case_id alt_count samples expected; do
  [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
  actual="${tmp_dir}/${case_id}.actual.tsv"
  "${cargo_run[@]}" annotate-core "${alt_count}" "${repo_root}/${samples}" >"${actual}"
  cmp -s "${repo_root}/${expected}" "${actual}" || {
    echo "[hc-full-parity-i1-core] mismatch ${case_id}" >&2
    diff -u "${repo_root}/${expected}" "${actual}" >&2 || true
    exit 1
  }
  echo "[hc-full-parity-i1-core] ok ${case_id}"
done <"${cases_tsv}"

echo "[hc-full-parity-i1-core] all cases match goldens."
