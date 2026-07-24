#!/usr/bin/env bash
# Phase F.2 gate: GATK native Log10PairHMM (Java) vs frozen goldens.
# Rust `pairhmm-native-likelihoods` is reserved for a future native port; L2 uses the same Java oracle.
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=lib_pinned_gatk.sh
source "${script_dir}/lib_pinned_gatk.sh"

repo_root="${GATK_RS_REPO_ROOT}"
cd "${repo_root}"

if [[ "${PARITY_SKIP_HC_FULL_F2:-0}" == "1" ]]; then
  echo "[hc-full-parity-f2] skipped (PARITY_SKIP_HC_FULL_F2=1)"
  exit 0
fi

cases_tsv="${repo_root}/parity/fixtures/hc-full-parity/f2/cases.tsv"
report_dir="${repo_root}/parity/reports"
tmp_dir="${report_dir}/hc-full-parity-f2-tmp"
mkdir -p "${tmp_dir}"

dump_java="${script_dir}/run_hc_full_parity_java_dump.sh"

"${script_dir}/run_hc_full_parity_java_compile.sh" >/dev/null

while IFS=$'\t' read -r case_id cases expected; do
  [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
  actual="${tmp_dir}/${case_id}.actual.tsv"
  "${dump_java}" pairhmm-native-likelihoods "${repo_root}/${cases}" >"${actual}" 2>"${tmp_dir}/${case_id}.stderr"
  if ! cmp -s "${repo_root}/${expected}" "${actual}"; then
    echo "[hc-full-parity-f2] mismatch case=${case_id}" >&2
    diff -u "${repo_root}/${expected}" "${actual}" >&2 || true
    exit 1
  fi
  echo "[hc-full-parity-f2] ok ${case_id}"
done <"${cases_tsv}"

echo "[hc-full-parity-f2] all cases match Java native goldens."
