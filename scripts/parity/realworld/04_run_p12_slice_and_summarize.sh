#!/usr/bin/env bash
# Real-world playbook — step 04: full Java + Rust HC on NA12878 slice + P12 JSON (long: Docker + compile).
# Prerequisites: 02 + 03 complete; set P12_REFERENCE or use default from 03 output.
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "${repo_root}"

if [[ -z "${P12_INTERVAL:-}" ]]; then
  echo "ERROR: set a bounded HC interval, e.g. export P12_INTERVAL='20:413419-463418'" >&2
  echo "  (Without -L, Java/Rust may attempt the whole BAM and run for a very long time.)" >&2
  exit 2
fi

ref="${P12_REFERENCE:-${repo_root}/parity/realworld/assets/hs37d5.simple.fa}"
export P12_REFERENCE="${ref}"
echo "=== 04_run_p12_slice_and_summarize P12_REFERENCE=${P12_REFERENCE} P12_INTERVAL=${P12_INTERVAL} ==="

./scripts/parity/run_p12_realworld_na12878_20k.sh
echo "04: reports under parity/reports/p12_realworld_na12878_20k.*"
