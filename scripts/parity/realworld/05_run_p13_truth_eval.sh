#!/usr/bin/env bash
# Real-world playbook — step 05: P13 metrics vs GIAB (needs 03 + 04 artifacts).
# Exports eval interval from P12_INTERVAL when set; must match the HC window used in step 04.
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "${repo_root}"

export P13_TRUTH_VCF="${P13_TRUTH_VCF:-${repo_root}/parity/realworld/assets/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz}"
export P13_REGIONS_BED="${P13_REGIONS_BED:-${repo_root}/parity/realworld/assets/HG001_GRCh37_1_22_v4.2.1_benchmark.bed}"
export P13_CHROM="${P13_CHROM:-20}"
# Prefer explicit P13_EVAL_INTERVAL, else reuse the same window as P12/P13 full harness
export P13_EVAL_INTERVAL="${P13_EVAL_INTERVAL:-${P12_INTERVAL:-}}"

echo "=== 05_run_p13_truth_eval ==="
echo "P13_TRUTH_VCF=${P13_TRUTH_VCF}"
echo "P13_EVAL_INTERVAL=${P13_EVAL_INTERVAL:-<empty = whole chrom in BED scope>}"
./scripts/parity/run_p13_realworld_truth_eval.sh
echo "05: see parity/reports/p13_realworld_truth_eval.json"
