#!/usr/bin/env bash
# CombineGVCFs parity: Java GATK 4.4 vs gatk-rs on a tiny synthetic mini cohort.
#
# Usage (Mac-local):
#   ./scripts/parity/run_combine_gvcfs_parity.sh
#
# Env:
#   COMBINE_SKIP_JAVA=1  — reuse existing Java output
#   COMBINE_SKIP_RUST=1  — reuse existing Rust output
#   GATK_JAR / gatk / GATK_DOCKER_IMAGE — via lib_pinned_gatk.sh + run_java_gatk.sh
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
log="${report_dir}/combine_gvcfs_${stamp}.log"
java_out="${report_dir}/combine_gvcfs.java.g.vcf"
rust_out="${report_dir}/combine_gvcfs.rust.g.vcf"
exec > >(tee -a "${log}") 2>&1

echo "=== CombineGVCFs parity ${stamp} ==="
echo "ref=${ref}"
echo "inputs: ${s1} ${s2}"

if [[ ! -f "${ref}.fai" ]]; then
  if command -v samtools >/dev/null 2>&1; then
    samtools faidx "${ref}"
  else
    echo "[combine-gvcfs] warning: no samtools; Java may require ${ref}.fai" >&2
  fi
fi
if [[ ! -f "${fixture_dir}/ref.dict" ]]; then
  echo "[combine-gvcfs] missing ${fixture_dir}/ref.dict (required by Java GATK)" >&2
  exit 1
fi

set +e
if [[ "${COMBINE_SKIP_JAVA:-0}" == "1" && -s "${java_out}" ]]; then
  echo "[combine-gvcfs] skipping Java (COMBINE_SKIP_JAVA=1)"
  java_exit=0
else
  echo "[combine-gvcfs] indexing input gVCFs for Java (IndexFeatureFile)…"
  "${repo_root}/scripts/parity/run_java_gatk.sh" \
    "${report_dir}/combine_gvcfs.index1.stdout" \
    IndexFeatureFile -I "${s1}" || true
  "${repo_root}/scripts/parity/run_java_gatk.sh" \
    "${report_dir}/combine_gvcfs.index2.stdout" \
    IndexFeatureFile -I "${s2}" || true

  echo "[combine-gvcfs] running Java CombineGVCFs…"
  : >"${report_dir}/combine_gvcfs.java.stdout"
  "${repo_root}/scripts/parity/run_java_gatk.sh" \
    "${report_dir}/combine_gvcfs.java.stdout" \
    CombineGVCFs \
    -R "${ref}" \
    -V "${s1}" \
    -V "${s2}" \
    -O "${java_out}"
  java_exit=$?
fi

target_dir="${CARGO_TARGET_DIR:-${repo_root}/target}"
rust_bin="${target_dir}/debug/gatk-rs"
if [[ ! -x "${rust_bin}" ]]; then
  rust_bin="${target_dir}/release/gatk-rs"
fi
if [[ ! -x "${rust_bin}" ]]; then
  echo "[combine-gvcfs] building gatk-rs (debug)…"
  export CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-1}"
  cargo build -p gatk-cli --bin gatk-rs
  rust_bin="${target_dir}/debug/gatk-rs"
fi

if [[ "${COMBINE_SKIP_RUST:-0}" == "1" && -s "${rust_out}" ]]; then
  echo "[combine-gvcfs] skipping Rust (COMBINE_SKIP_RUST=1)"
  rust_exit=0
else
  echo "[combine-gvcfs] running Rust CombineGVCFs…"
  "${rust_bin}" combine-gvcfs \
    -R "${ref}" \
    -V "${s1}" \
    -V "${s2}" \
    -O "${rust_out}"
  rust_exit=$?
fi
set -e

echo "[combine-gvcfs] java_exit=${java_exit} rust_exit=${rust_exit}"
if [[ "${java_exit}" -ne 0 ]]; then
  echo "[combine-gvcfs] Java failed; see ${report_dir}/combine_gvcfs.java.stdout" >&2
  exit "${java_exit}"
fi
if [[ "${rust_exit}" -ne 0 ]]; then
  echo "[combine-gvcfs] Rust failed" >&2
  exit "${rust_exit}"
fi

python3 "${repo_root}/scripts/parity/compare_combine_gvcfs.py" \
  --java "${java_out}" \
  --rust "${rust_out}" \
  --label combine-gvcfs-mini
cmp_exit=$?
echo "[combine-gvcfs] compare_exit=${cmp_exit}"
echo "[combine-gvcfs] log=${log}"
exit "${cmp_exit}"
