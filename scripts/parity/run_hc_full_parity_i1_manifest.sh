#!/usr/bin/env bash
# Phase I.0.1 — parity v1 annotation manifest dump.
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

expected="${repo_root}/parity/fixtures/hc-full-parity/i1/expected/manifest.tsv"
tmp_dir="${repo_root}/parity/reports/hc-full-parity-i1-tmp"
mkdir -p "${tmp_dir}"
actual="${tmp_dir}/manifest.actual.tsv"

cargo_run=(cargo run -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)
[[ "${PARITY_RUST_PROFILE:-dev}" == "release" ]] && cargo_run=(cargo run --release -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)

"${cargo_run[@]}" annotation-manifest >"${actual}"
cmp -s "${expected}" "${actual}" || {
  echo "[hc-full-parity-i1-manifest] mismatch" >&2
  diff -u "${expected}" "${actual}" >&2 || true
  exit 1
}
echo "[hc-full-parity-i1-manifest] ok"
