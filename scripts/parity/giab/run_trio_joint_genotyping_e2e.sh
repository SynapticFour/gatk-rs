#!/usr/bin/env bash
# End-to-end joint genotyping proof:
#   per-sample HaplotypeCaller (-ERC GVCF) → CombineGVCFs → GenotypeGVCFs
# for GIAB Ashkenazi trio HG002/HG003/HG004 on a small interval, scored with
# hap.py (via gatk-rs-equiv when available) against trio truth AND compared to
# the identical Java GATK 4.4 path.
#
# Usage:
#   # Smoke (synthetic mini cohort — no BAM download; always runnable):
#   TRIO_E2E_MODE=smoke ./scripts/parity/giab/run_trio_joint_genotyping_e2e.sh
#
#   # Real GIAB trio (requires staged BAMs + truth):
#   TRIO_E2E_MODE=giab \
#     TRIO_HG002_BAM=... TRIO_HG003_BAM=... TRIO_HG004_BAM=... \
#     TRIO_REFERENCE=... TRIO_INTERVAL=20:1000000-1050000 \
#     TRIO_TRUTH_VCF=... TRIO_TRUTH_BED=... \
#     ./scripts/parity/giab/run_trio_joint_genotyping_e2e.sh
#
# Env knobs:
#   TRIO_E2E_MODE=smoke|giab   (default: smoke)
#   TRIO_SKIP_JAVA=1 / TRIO_SKIP_RUST=1
#   TRIO_STAND_CALL_CONF=30
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "${repo_root}"
# shellcheck source=../lib_pinned_gatk.sh
source "${repo_root}/scripts/parity/lib_pinned_gatk.sh"

mode="${TRIO_E2E_MODE:-smoke}"
report_root="${repo_root}/parity/reports/trio_joint_e2e"
stamp="$(date -u +%Y%m%dT%H%M%SZ)"
run_dir="${report_root}/${stamp}_${mode}"
mkdir -p "${run_dir}"
log="${run_dir}/run.log"
exec > >(tee -a "${log}") 2>&1

echo "=== Trio joint genotyping E2E (${mode}) ${stamp} ==="

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

stand_call_conf="${TRIO_STAND_CALL_CONF:-30}"

run_smoke() {
  local fixture="${repo_root}/parity/combine_gvcfs/mini"
  local ref="${fixture}/ref.fa"
  local s1="${fixture}/sample1.g.vcf"
  local s2="${fixture}/sample2.g.vcf"
  local java_c="${run_dir}/java.combined.g.vcf"
  local rust_c="${run_dir}/rust.combined.g.vcf"
  local java_v="${run_dir}/java.final.vcf"
  local rust_v="${run_dir}/rust.final.vcf"

  "${repo_root}/scripts/parity/run_java_gatk.sh" \
    "${run_dir}/index1.stdout" IndexFeatureFile -I "${s1}" || true
  "${repo_root}/scripts/parity/run_java_gatk.sh" \
    "${run_dir}/index2.stdout" IndexFeatureFile -I "${s2}" || true

  if [[ "${TRIO_SKIP_JAVA:-0}" != "1" ]]; then
    "${repo_root}/scripts/parity/run_java_gatk.sh" \
      "${run_dir}/java.combine.stdout" \
      CombineGVCFs -R "${ref}" -V "${s1}" -V "${s2}" -O "${java_c}"
    "${repo_root}/scripts/parity/run_java_gatk.sh" \
      "${run_dir}/java.gg.stdout" \
      GenotypeGVCFs -R "${ref}" -V "${java_c}" -O "${java_v}" \
      --standard-min-confidence-threshold-for-calling "${stand_call_conf}"
  fi
  if [[ "${TRIO_SKIP_RUST:-0}" != "1" ]]; then
    "${rust_bin}" combine-gvcfs -R "${ref}" -V "${s1}" -V "${s2}" -O "${rust_c}"
    "${rust_bin}" genotype-gvcfs -R "${ref}" -V "${rust_c}" -O "${rust_v}" \
      --stand-call-conf "${stand_call_conf}"
  fi

  python3 "${repo_root}/scripts/parity/compare_genotype_gvcfs.py" \
    --java "${java_v}" \
    --rust "${rust_v}" \
    --label trio-e2e-smoke
}

run_giab() {
  local ref="${TRIO_REFERENCE:?set TRIO_REFERENCE}"
  local interval="${TRIO_INTERVAL:?set TRIO_INTERVAL (e.g. 20:1000000-1050000)}"
  local bam2="${TRIO_HG002_BAM:?set TRIO_HG002_BAM}"
  local bam3="${TRIO_HG003_BAM:?set TRIO_HG003_BAM}"
  local bam4="${TRIO_HG004_BAM:?set TRIO_HG004_BAM}"
  local truth_vcf="${TRIO_TRUTH_VCF:-}"
  local truth_bed="${TRIO_TRUTH_BED:-}"

  for b in "${bam2}" "${bam3}" "${bam4}" "${ref}"; do
    [[ -f "${b}" ]] || { echo "missing input: ${b}" >&2; exit 1; }
  done

  local samples=(HG002 HG003 HG004)
  local bams=("${bam2}" "${bam3}" "${bam4}")

  # --- Rust path ---
  local rust_gvcfs=()
  if [[ "${TRIO_SKIP_RUST:-0}" != "1" ]]; then
    for i in 0 1 2; do
      local g="${run_dir}/rust.${samples[$i]}.g.vcf"
      echo "[trio-e2e] Rust HC GVCF ${samples[$i]}…"
      "${rust_bin}" haplotypecaller \
        -R "${ref}" -I "${bams[$i]}" -O "${g}" -L "${interval}" \
        --emit-ref-confidence GVCF
      rust_gvcfs+=("${g}")
    done
    local rust_c="${run_dir}/rust.combined.g.vcf"
    local rust_v="${run_dir}/rust.final.vcf"
    "${rust_bin}" combine-gvcfs -R "${ref}" \
      -V "${rust_gvcfs[0]}" -V "${rust_gvcfs[1]}" -V "${rust_gvcfs[2]}" \
      -O "${rust_c}"
    "${rust_bin}" genotype-gvcfs -R "${ref}" -V "${rust_c}" -O "${rust_v}" \
      --stand-call-conf "${stand_call_conf}"
  fi

  # --- Java path ---
  local java_gvcfs=()
  if [[ "${TRIO_SKIP_JAVA:-0}" != "1" ]]; then
    for i in 0 1 2; do
      local g="${run_dir}/java.${samples[$i]}.g.vcf"
      echo "[trio-e2e] Java HC GVCF ${samples[$i]}…"
      "${repo_root}/scripts/parity/run_java_gatk.sh" \
        "${run_dir}/java.hc.${samples[$i]}.stdout" \
        HaplotypeCaller -R "${ref}" -I "${bams[$i]}" -O "${g}" -L "${interval}" \
        -ERC GVCF --verbosity ERROR
      java_gvcfs+=("${g}")
      "${repo_root}/scripts/parity/run_java_gatk.sh" \
        "${run_dir}/java.idx.${samples[$i]}.stdout" \
        IndexFeatureFile -I "${g}" || true
    done
    local java_c="${run_dir}/java.combined.g.vcf"
    local java_v="${run_dir}/java.final.vcf"
    "${repo_root}/scripts/parity/run_java_gatk.sh" \
      "${run_dir}/java.combine.stdout" \
      CombineGVCFs -R "${ref}" \
      -V "${java_gvcfs[0]}" -V "${java_gvcfs[1]}" -V "${java_gvcfs[2]}" \
      -O "${java_c}"
    "${repo_root}/scripts/parity/run_java_gatk.sh" \
      "${run_dir}/java.gg.stdout" \
      GenotypeGVCFs -R "${ref}" -V "${java_c}" -O "${java_v}" \
      --standard-min-confidence-threshold-for-calling "${stand_call_conf}"
  fi

  echo "[trio-e2e] comparing Java vs Rust final VCFs…"
  python3 "${repo_root}/scripts/parity/compare_genotype_gvcfs.py" \
    --java "${run_dir}/java.final.vcf" \
    --rust "${run_dir}/rust.final.vcf" \
    --label trio-e2e-giab-java-rust \
    --qual-tol 50

  if [[ -n "${truth_vcf}" && -f "${truth_vcf}" ]]; then
    local happy_bin="${HAPPY_BIN:-}"
    if [[ -z "${happy_bin}" ]] && command -v hap.py >/dev/null 2>&1; then
      happy_bin="hap.py"
    fi
    if [[ -n "${happy_bin}" ]]; then
      echo "[trio-e2e] hap.py vs truth (Rust final VCF)…"
      mkdir -p "${run_dir}/happy_rust" "${run_dir}/happy_java"
      local bed_args=()
      if [[ -n "${truth_bed}" && -f "${truth_bed}" ]]; then
        bed_args+=(-f "${truth_bed}")
      fi
      "${happy_bin}" "${truth_vcf}" "${run_dir}/rust.final.vcf" \
        -r "${ref}" "${bed_args[@]}" -o "${run_dir}/happy_rust/prefix" \
        --threads "${TRIO_HAPPY_THREADS:-2}" || {
          echo "[trio-e2e] warning: hap.py (Rust) failed"
        }
      echo "[trio-e2e] hap.py vs truth (Java final VCF)…"
      "${happy_bin}" "${truth_vcf}" "${run_dir}/java.final.vcf" \
        -r "${ref}" "${bed_args[@]}" -o "${run_dir}/happy_java/prefix" \
        --threads "${TRIO_HAPPY_THREADS:-2}" || {
          echo "[trio-e2e] warning: hap.py (Java) failed"
        }
      # Also stash a pointer for tools/equivalence consumers.
      {
        echo "rust_happy=${run_dir}/happy_rust"
        echo "java_happy=${run_dir}/happy_java"
        echo "rust_vcf=${run_dir}/rust.final.vcf"
        echo "java_vcf=${run_dir}/java.final.vcf"
      } >"${run_dir}/equiv_paths.txt"
    else
      echo "[trio-e2e] hap.py not found (set HAPPY_BIN); Java↔Rust compare only"
    fi
  else
    echo "[trio-e2e] TRIO_TRUTH_VCF unset/missing — Java↔Rust compare only"
  fi
}

case "${mode}" in
  smoke) run_smoke ;;
  giab) run_giab ;;
  *)
    echo "unknown TRIO_E2E_MODE=${mode} (use smoke|giab)" >&2
    exit 2
    ;;
esac

echo "[trio-e2e] done run_dir=${run_dir}"
echo "${run_dir}" >"${report_root}/latest.txt"
