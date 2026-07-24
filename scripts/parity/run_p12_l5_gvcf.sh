#!/usr/bin/env bash
# P12 L5.2 — active gVCF on NA12878 chr2 slice vs Java (-ERC GVCF).
#
# Usage:
#   export P12_REFERENCE="$PWD/parity/realworld/assets/hs37d5.simple.fa"
#   ./scripts/parity/run_p12_l5_gvcf.sh
#
# Logs: parity/reports/p12_l5_gvcf_<timestamp>.log
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"

if [[ -z "${P12_REFERENCE:-}" ]]; then
  export P12_REFERENCE="${repo_root}/parity/realworld/assets/hs37d5.simple.fa"
fi
if [[ ! -f "${P12_REFERENCE}" ]]; then
  echo "P12_REFERENCE not found: ${P12_REFERENCE}" >&2
  exit 1
fi

bam="${P12_BAM:-${repo_root}/parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam}"
if [[ ! -f "${bam}" ]]; then
  echo "P12 BAM not found: ${bam}" >&2
  echo "Stage with: ./scripts/parity/realworld/02_stage_na12878_20k_bam.sh" >&2
  exit 1
fi

interval="${P12_INTERVAL:-2:92300000-92350000}"
report_dir="${repo_root}/parity/reports"
mkdir -p "${report_dir}"

stamp="$(date -u +%Y%m%dT%H%M%SZ)"
log="${report_dir}/p12_l5_gvcf_${stamp}.log"
canonical="${report_dir}/p12_l5_gvcf_canonical.log"
exec > >(tee -a "${log}") 2>&1

java_gvcf="${report_dir}/p12_l5_gvcf.java.g.vcf"
rust_gvcf="${report_dir}/p12_l5_gvcf.rust.g.vcf"
json_out="${report_dir}/p12_l5_gvcf.json"
md_out="${report_dir}/p12_l5_gvcf.md"

target_dir="${CARGO_TARGET_DIR:-${repo_root}/target}"
export RAYON_NUM_THREADS="${RAYON_NUM_THREADS:-2}"
export CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-1}"

echo "=== P12 L5.2 gVCF ${stamp} ==="
echo "[p12-l5-gvcf] tip: P12_SKIP_JAVA=1 P12_SKIP_RUST=1 skips HC (~1s compare). Full Rust re-run ~15–20 min."
echo "P12_REFERENCE=${P12_REFERENCE}"
echo "P12_BAM=${bam}"
echo "P12_INTERVAL=${interval}"

set +e
if [[ "${P12_SKIP_JAVA:-0}" == "1" && -s "${java_gvcf}" ]]; then
  echo "[p12-l5-gvcf] skipping Java HC (P12_SKIP_JAVA=1, existing ${java_gvcf})"
  java_exit=0
else
  echo "[p12-l5-gvcf] running Java HaplotypeCaller -ERC GVCF via Docker (often 5–15 min on arm64)…"
  docker run --rm --platform "${GATK_DOCKER_PLATFORM:-linux/amd64}" \
    -v "${repo_root}:${repo_root}" \
    -w "${repo_root}" \
    "${GATK_DOCKER_IMAGE:-us.gcr.io/broad-gatk/gatk:4.4.0.0}" \
    gatk HaplotypeCaller \
    -R "${P12_REFERENCE}" \
    -I "${bam}" \
    -O "${java_gvcf}" \
    -L "${interval}" \
    -ERC GVCF \
    --verbosity ERROR
  java_exit=$?
fi

rust_bin="${target_dir}/release/gatk-rs"
if [[ ! -x "${rust_bin}" ]]; then
  echo "[p12-l5-gvcf] building release gatk-rs…"
  ./scripts/parity/build_gatk_rs_release.sh
fi

if [[ "${P12_SKIP_RUST:-0}" == "1" && -s "${rust_gvcf}" ]]; then
  echo "[p12-l5-gvcf] skipping Rust HC (P12_SKIP_RUST=1, existing ${rust_gvcf})"
  rust_exit=0
else
  echo "[p12-l5-gvcf] running Rust HaplotypeCaller --output-mode GVCF (~2–5 min)…"
  "${rust_bin}" HaplotypeCaller \
    -R "${P12_REFERENCE}" \
    -I "${bam}" \
    -O "${rust_gvcf}" \
    -L "${interval}" \
    --output-mode GVCF
  rust_exit=$?
fi
set -e

if [[ "${java_exit}" -ne 0 || "${rust_exit}" -ne 0 ]]; then
  echo "[p12-l5-gvcf] tool error: java_exit=${java_exit} rust_exit=${rust_exit}" >&2
  exit 1
fi

python3 "${repo_root}/scripts/parity/compare_p12_l5_gvcf.py" \
  --java "${java_gvcf}" \
  --rust "${rust_gvcf}" \
  --block-contract "${P12_BLOCK_CONTRACT:-semantic}" \
  --json-out "${json_out}" \
  --md-out "${md_out}"
compare_exit=$?

cp -f "${log}" "${canonical}"

echo ""
if [[ "${compare_exit}" -eq 0 ]]; then
  echo "=== P12 L5.2 gVCF PASS (variant + block) ==="
else
  status="$(python3 -c "import json; print(json.load(open('${json_out}'))['status'])")"
  if [[ "${status}" == "variant_pass" ]]; then
    echo "=== P12 L5.2 gVCF VARIANT PASS (block boundaries open) ==="
  else
    echo "=== P12 L5.2 gVCF FAIL ==="
    exit 1
  fi
fi
echo "log: ${log}"
echo "canonical: ${canonical}"
echo "json: ${json_out}"
echo "md: ${md_out}"
