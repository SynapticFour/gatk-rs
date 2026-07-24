#!/usr/bin/env bash
# L5.4 — chr2 2k VCF scale gate (subset of 20k slice).
#
# Usage:
#   export P12_REFERENCE="$PWD/parity/realworld/assets/hs37d5.simple.fa"
#   ./scripts/parity/run_p12_l5_chr2_scale.sh
#
# Compare-only: P12_SKIP_JAVA=1 P12_SKIP_RUST=1 ./scripts/parity/run_p12_l5_chr2_scale.sh
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"

if [[ -z "${P12_REFERENCE:-}" ]]; then
  export P12_REFERENCE="${repo_root}/parity/realworld/assets/hs37d5.simple.fa"
fi

bam="${P12_BAM:-${repo_root}/parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam}"
interval="${P12_SCALE_INTERVAL:-2:92300000-92302000}"
report_dir="${repo_root}/parity/reports"
mkdir -p "${report_dir}"

java_vcf="${report_dir}/p12_l5_chr2_scale.java.vcf"
rust_vcf="${report_dir}/p12_l5_chr2_scale.rust.vcf"
json_out="${report_dir}/p12_l5_chr2_scale.json"
md_out="${report_dir}/p12_l5_chr2_scale.md"

target_dir="${CARGO_TARGET_DIR:-${repo_root}/target}"
export RAYON_NUM_THREADS="${RAYON_NUM_THREADS:-2}"

echo "=== P12 L5.4 chr2 scale ${interval} ==="
echo "P12_REFERENCE=${P12_REFERENCE}"
echo "P12_BAM=${bam}"

set +e
if [[ "${P12_SKIP_JAVA:-0}" == "1" && -s "${java_vcf}" ]]; then
  echo "[p12-l5-scale] skipping Java HC"
  java_exit=0
else
  docker run --rm --platform "${GATK_DOCKER_PLATFORM:-linux/amd64}" \
    -v "${repo_root}:${repo_root}" \
    -w "${repo_root}" \
    "${GATK_DOCKER_IMAGE:-us.gcr.io/broad-gatk/gatk:4.4.0.0}" \
    gatk HaplotypeCaller \
    -R "${P12_REFERENCE}" \
    -I "${bam}" \
    -O "${java_vcf}" \
    -L "${interval}" \
    -ERC GVCF \
    --verbosity ERROR
  java_exit=$?
fi

rust_bin="${target_dir}/release/gatk-rs"
if [[ ! -x "${rust_bin}" ]]; then
  ./scripts/parity/build_gatk_rs_release.sh
fi

if [[ "${P12_SKIP_RUST:-0}" == "1" && -s "${rust_vcf}" ]]; then
  echo "[p12-l5-scale] skipping Rust HC"
  rust_exit=0
else
  "${rust_bin}" HaplotypeCaller \
    -R "${P12_REFERENCE}" \
    -I "${bam}" \
    -O "${rust_vcf}" \
    -L "${interval}" \
    --output-mode GVCF
  rust_exit=$?
fi
set -e

if [[ "${java_exit}" -ne 0 || "${rust_exit}" -ne 0 ]]; then
  echo "[p12-l5-scale] tool error: java_exit=${java_exit} rust_exit=${rust_exit}" >&2
  exit 1
fi

python3 "${repo_root}/scripts/parity/compare_p12_l5_chr2_scale.py" \
  --java "${java_vcf}" \
  --rust "${rust_vcf}" \
  --json-out "${json_out}" \
  --md-out "${md_out}"
