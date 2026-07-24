#!/usr/bin/env bash
# GenotypeGVCFs parity: Java GATK 4.4 vs gatk-rs on CombineGVCFs mini cohort output.
#
# Usage:
#   ./scripts/parity/run_genotype_gvcfs_parity.sh
#
# Env:
#   GG_SKIP_JAVA=1 / GG_SKIP_RUST=1 — reuse existing outputs
#   GG_STAND_CALL_CONF — default 30
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
# shellcheck source=lib_pinned_gatk.sh
source "${repo_root}/scripts/parity/lib_pinned_gatk.sh"

fixture_dir="${repo_root}/parity/combine_gvcfs/mini"
ref="${fixture_dir}/ref.fa"
s1="${fixture_dir}/sample1.g.vcf"
s2="${fixture_dir}/sample2.g.vcf"
report_dir="${repo_root}/parity/reports"
mkdir -p "${report_dir}"

stamp="$(date -u +%Y%m%dT%H%M%SZ)"
log="${report_dir}/genotype_gvcfs_${stamp}.log"
java_combined="${report_dir}/gg_parity.java.combined.g.vcf"
rust_combined="${report_dir}/gg_parity.rust.combined.g.vcf"
java_out="${report_dir}/genotype_gvcfs.java.vcf"
rust_out="${report_dir}/genotype_gvcfs.rust.vcf"
stand_call_conf="${GG_STAND_CALL_CONF:-30}"
exec > >(tee -a "${log}") 2>&1

echo "=== GenotypeGVCFs parity ${stamp} ==="

target_dir="${CARGO_TARGET_DIR:-${repo_root}/target}"
rust_bin="${target_dir}/debug/gatk-rs"
if [[ ! -x "${rust_bin}" ]]; then
  rust_bin="${target_dir}/release/gatk-rs"
fi
if [[ ! -x "${rust_bin}" ]]; then
  export CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-1}"
  cargo build -p gatk-cli --bin gatk-rs
  rust_bin="${target_dir}/debug/gatk-rs"
fi

# Ensure Java indexes for Combine inputs.
"${repo_root}/scripts/parity/run_java_gatk.sh" \
  "${report_dir}/gg_parity.index1.stdout" \
  IndexFeatureFile -I "${s1}" || true
"${repo_root}/scripts/parity/run_java_gatk.sh" \
  "${report_dir}/gg_parity.index2.stdout" \
  IndexFeatureFile -I "${s2}" || true

set +e
if [[ "${GG_SKIP_JAVA:-0}" == "1" && -s "${java_out}" ]]; then
  echo "[gg-parity] skipping Java"
  java_exit=0
else
  echo "[gg-parity] Java CombineGVCFs…"
  "${repo_root}/scripts/parity/run_java_gatk.sh" \
    "${report_dir}/gg_parity.java.combine.stdout" \
    CombineGVCFs -R "${ref}" -V "${s1}" -V "${s2}" -O "${java_combined}"
  echo "[gg-parity] Java GenotypeGVCFs…"
  "${repo_root}/scripts/parity/run_java_gatk.sh" \
    "${report_dir}/gg_parity.java.gg.stdout" \
    GenotypeGVCFs -R "${ref}" -V "${java_combined}" -O "${java_out}" \
    --standard-min-confidence-threshold-for-calling "${stand_call_conf}"
  java_exit=$?
fi

if [[ "${GG_SKIP_RUST:-0}" == "1" && -s "${rust_out}" ]]; then
  echo "[gg-parity] skipping Rust"
  rust_exit=0
else
  echo "[gg-parity] Rust CombineGVCFs…"
  "${rust_bin}" combine-gvcfs -R "${ref}" -V "${s1}" -V "${s2}" -O "${rust_combined}"
  echo "[gg-parity] Rust GenotypeGVCFs…"
  "${rust_bin}" genotype-gvcfs -R "${ref}" -V "${rust_combined}" -O "${rust_out}" \
    --stand-call-conf "${stand_call_conf}"
  rust_exit=$?
fi
set -e

echo "[gg-parity] java_exit=${java_exit} rust_exit=${rust_exit}"
if [[ "${java_exit}" -ne 0 || "${rust_exit}" -ne 0 ]]; then
  exit 1
fi

python3 "${repo_root}/scripts/parity/compare_genotype_gvcfs.py" \
  --java "${java_out}" \
  --rust "${rust_out}" \
  --label genotype-gvcfs-mini
cmp_exit=$?
echo "[gg-parity] compare_exit=${cmp_exit} log=${log}"
exit "${cmp_exit}"
