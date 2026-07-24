#!/usr/bin/env bash
# Real-world playbook — step 02: download NA12878_20k b37 BAM/BAI only (no HC, no reference required).
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "${repo_root}"

echo "=== 02_stage_na12878_20k_bam ==="
export P12_REFERENCE=""
./scripts/parity/run_p12_realworld_na12878_20k.sh
echo "02_stage_na12878_20k_bam: done (see parity/realworld/na12878_20k_b37/)"
