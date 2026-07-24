#!/usr/bin/env bash
# Foundation evidence (01+08–11) + P13 refresh on current P12 VCFs + evidence markdown/json.
# Always refreshes P13 + report even if parity smoke fails (e.g. Java GATK not on PATH — use Docker per run_java_gatk.sh).
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../../.." && pwd)"
rd="${repo_root}/scripts/parity/realworld"
cd "${repo_root}"

export P12_INTERVAL="${P12_INTERVAL:-20:413419-463418}"

set +e
"${rd}/run_foundation_evidence.sh"
foundation_ec=$?
set -e

export P13_TRUTH_VCF="${P13_TRUTH_VCF:-${repo_root}/parity/realworld/assets/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz}"
export P13_REGIONS_BED="${P13_REGIONS_BED:-${repo_root}/parity/realworld/assets/HG001_GRCh37_1_22_v4.2.1_benchmark.bed}"
export P13_CHROM="${P13_CHROM:-20}"
export P13_EVAL_INTERVAL="${P13_EVAL_INTERVAL:-${P12_INTERVAL}}"

"${rd}/05_run_p13_truth_eval.sh"

python3 "${rd}/generate_evidence_report.py" "${repo_root}/parity/reports/realworld_run_manifest.json"

echo "run_full_evidence: see parity/reports/realworld_parity_evidence.md (foundation exit=${foundation_ec})"
exit "${foundation_ec}"
