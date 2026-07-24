#!/usr/bin/env bash
# J-D05 — DRAGEN mode branch scaffold (default off).
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

expected="${repo_root}/parity/fixtures/hc-full-parity/j-dragen/expected/default_off.tsv"
tmp_dir="${repo_root}/parity/reports/hc-full-parity-j-dragen-tmp"
mkdir -p "${tmp_dir}"
actual="${tmp_dir}/default_off.actual.tsv"

cargo_run=(cargo run -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)
[[ "${PARITY_RUST_PROFILE:-dev}" == "release" ]] && cargo_run=(cargo run --release -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump --)

"${cargo_run[@]}" dragen-mode-branch >"${actual}" 2>"${tmp_dir}/default_off.stderr"
cmp -s "${expected}" "${actual}" || {
  echo "[hc-full-parity-j-dragen] mismatch" >&2
  diff -u "${expected}" "${actual}" >&2 || true
  exit 1
}
echo "[hc-full-parity-j-dragen] ok default_off"
