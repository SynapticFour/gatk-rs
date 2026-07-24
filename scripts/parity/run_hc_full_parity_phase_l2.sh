#!/usr/bin/env bash
# Phase A P1 bundle: Java L2 dumps present + Rust vs Java comparison (strict by default; PARITY_HC_FULL_L2_STRICT=0 to allow mismatches).
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"
cd "${repo_root}"

if [[ "${PARITY_SKIP_HC_FULL_L2:-0}" == "1" ]]; then
  echo "[hc-full-parity-phase-l2] skipped"
  exit 0
fi

pin_short="$(grep GATK_PINNED_SHA_SHORT docs/GATK_PINNED.env | cut -d= -f2)"
java_probe="${repo_root}/parity/fixtures/hc-full-parity/java_dumps/b1/chr1_5_15_default_${pin_short}.tsv"
if [[ ! -f "${java_probe}" ]]; then
  echo "[hc-full-parity-phase-l2] missing Java dumps; run ./scripts/parity/run_hc_full_parity_java_refresh.sh" >&2
  exit 2
fi

"${script_dir}/run_hc_full_parity_l2.sh"
