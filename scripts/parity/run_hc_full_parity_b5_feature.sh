#!/usr/bin/env bash
# Phase B.5.4 — FeatureContext per assembly region (padded span).
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

if [[ "${PARITY_SKIP_HC_FULL_B5_FEATURE:-0}" == "1" ]]; then
  echo "[hc-full-parity-b5-feature] skipped (PARITY_SKIP_HC_FULL_B5_FEATURE=1)"
  exit 0
fi

cases_tsv="${repo_root}/parity/fixtures/hc-full-parity/b5-feature/cases.tsv"
report_dir="${repo_root}/parity/reports"
tmp_dir="${report_dir}/hc-full-parity-b5-feature-tmp"
mkdir -p "${tmp_dir}"

profile="${PARITY_RUST_PROFILE:-dev}"
cargo_run=(cargo run -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)
if [[ "${profile}" == "release" ]]; then
  cargo_run=(cargo run --release -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)
fi

while IFS=$'\t' read -r case_id ref bam interval padding features_vcf expected; do
  [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
  actual="${tmp_dir}/${case_id}.actual.tsv"
  set +e
  args=(assembly-region-features "${repo_root}/${ref}" "${repo_root}/${bam}" "${interval}")
  if [[ -n "${padding}" && "${padding}" != "-" ]]; then
    args+=("${padding}")
  elif [[ -n "${features_vcf}" && "${features_vcf}" != "-" ]]; then
    args+=("-")
  fi
  if [[ -n "${features_vcf}" && "${features_vcf}" != "-" ]]; then
    args+=("${repo_root}/${features_vcf}")
  fi
  "${cargo_run[@]}" "${args[@]}" >"${actual}" 2>"${tmp_dir}/${case_id}.stderr"
  dump_ec=$?
  set -e
  if [[ "${dump_ec}" -ne 0 ]]; then
    echo "[hc-full-parity-b5-feature] dump failed (case=${case_id})" >&2
    cat "${tmp_dir}/${case_id}.stderr" >&2
    exit 1
  fi
  if ! cmp -s "${repo_root}/${expected}" "${actual}"; then
    echo "[hc-full-parity-b5-feature] mismatch case=${case_id}" >&2
    diff -u "${repo_root}/${expected}" "${actual}" >&2 || true
    exit 1
  fi
  echo "[hc-full-parity-b5-feature] ok ${case_id}"
done <"${cases_tsv}"

echo "[hc-full-parity-b5-feature] all cases match goldens."
