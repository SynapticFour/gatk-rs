#!/usr/bin/env bash
# I-D01 — standard annotation plugins (FS, QD, BaseQRankSum).
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

cases_tsv="${repo_root}/parity/fixtures/hc-full-parity/i1-standard/cases.tsv"
tmp_dir="${repo_root}/parity/reports/hc-full-parity-i1-standard-tmp"
mkdir -p "${tmp_dir}"

cargo_run=(cargo run -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)
[[ "${PARITY_RUST_PROFILE:-dev}" == "release" ]] && cargo_run=(cargo run --release -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)

while IFS=$'\t' read -r case_id ref_fw ref_rv alt_fw alt_rv qual dp ref_bqs alt_bqs ref_pos alt_pos ref_mq alt_mq expected; do
  [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
  actual="${tmp_dir}/${case_id}.actual.tsv"
  "${cargo_run[@]}" standard-annotations "${ref_fw}" "${ref_rv}" "${alt_fw}" "${alt_rv}" "${qual}" "${dp}" "${ref_bqs}" "${alt_bqs}" "${ref_pos}" "${alt_pos}" "${ref_mq}" "${alt_mq}" >"${actual}" 2>"${tmp_dir}/${case_id}.stderr"
  cmp -s "${repo_root}/${expected}" "${actual}" || {
    echo "[hc-full-parity-i1-standard] mismatch case=${case_id}" >&2
    diff -u "${repo_root}/${expected}" "${actual}" >&2 || true
    exit 1
  }
  echo "[hc-full-parity-i1-standard] ok ${case_id}"
done <"${cases_tsv}"

echo "[hc-full-parity-i1-standard] all cases match goldens."
