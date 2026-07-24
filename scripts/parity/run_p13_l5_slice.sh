#!/usr/bin/env bash
# P7 — truth eval on P12 L5 slice using gVCF variant rows.
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"

report_dir="${repo_root}/parity/reports"
truth="${P13_TRUTH_VCF:-${repo_root}/parity/realworld/assets/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz}"
java_vcf="${P13_JAVA_VCF:-${report_dir}/p12_l5_gvcf.java.g.vcf}"
rust_vcf="${P13_RUST_VCF:-${report_dir}/p12_l5_gvcf.rust.g.vcf}"
eval_interval="${P13_EVAL_INTERVAL:-2:92300000-92350000}"

export P13_TRUTH_VCF="${truth}"
export P13_JAVA_VCF="${java_vcf}"
export P13_RUST_VCF="${rust_vcf}"
export P13_EVAL_INTERVAL="${eval_interval}"
export P13_REGIONS_BED="${P13_REGIONS_BED:-${repo_root}/parity/realworld/assets/HG001_GRCh37_1_22_v4.2.1_benchmark.bed}"
export P12_INTERVAL="${eval_interval}"

./scripts/parity/run_p13_realworld_truth_eval.sh
