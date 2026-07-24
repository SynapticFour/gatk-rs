#!/usr/bin/env bash
# Real-world playbook — step 06: re-run Rust HC only (reuse cached Java VCF from step 04).
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "${repo_root}"

if [[ -z "${P12_INTERVAL:-}" ]]; then
  echo "ERROR: set P12_INTERVAL to match the Java run window, e.g. export P12_INTERVAL='20:413419-463418'" >&2
  exit 2
fi

export P12_REFERENCE="${P12_REFERENCE:-${repo_root}/parity/realworld/assets/hs37d5.simple.fa}"
echo "=== 06_run_rust_only_resync P12_REFERENCE=${P12_REFERENCE} P12_INTERVAL=${P12_INTERVAL} ==="
./scripts/parity/run_p12_rust_only_na12878_20k.sh
echo "06: refreshed parity/reports/p12_realworld_na12878_20k.rust.vcf + json"
