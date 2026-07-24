#!/usr/bin/env bash
# Real-world playbook — step 03: hs37d5 (simple.fa + faidx + dict) + GIAB benchmark VCF/BED.
# Uses the same logic as run_p12_p13_realworld_full.sh but stops before interval/HC.
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "${repo_root}"

echo "=== 03_stage_reference_and_truth (REALWORLD_STOP_AFTER_ASSETS=1) ==="
export REALWORLD_STOP_AFTER_ASSETS=1
./scripts/parity/run_p12_p13_realworld_full.sh
echo "03_stage_reference_and_truth: done (parity/realworld/assets/)"
