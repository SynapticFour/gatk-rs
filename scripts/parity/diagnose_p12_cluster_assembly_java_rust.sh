#!/usr/bin/env bash
# P12 cluster: Java vs Rust assembly dumps on the same interval (ASM-1 side-by-side).
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"

REF="${P12_REFERENCE:-${repo_root}/parity/realworld/assets/hs37d5.simple.fa}"
BAM="${P12_BAM:-${repo_root}/parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam}"
# Rust walker: one active region 92307200–92307400. Java walker: active 92307228–92307400 only.
# Use Java-active sub-interval for aligned graph/haplotype dumps (see assembly-regions).
INT="${P12_CLUSTER_INTERVAL:-2:92307228-92307400}"
RUST_FULL_INTERVAL="${P12_RUST_FULL_INTERVAL:-2:92307200-92307400}"
OUT="${P12_SIDE_BY_SIDE_DIR:-${repo_root}/parity/reports/p12_cluster_side_by_side}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"
RUST_DUMP="${PARITY_DUMP_BIN:-${CARGO_TARGET_DIR}/release/examples/hc_full_parity_gate_dump}"
JAVA_DUMP="${repo_root}/scripts/parity/run_hc_full_parity_java_dump.sh"

mkdir -p "${OUT}"

if [[ ! -x "${RUST_DUMP}" ]]; then
  echo "Build Rust dump: CARGO_TARGET_DIR=${repo_root}/target cargo build --release -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump"
  exit 1
fi

echo "[side-by-side] aligned_interval=${INT} (Java-active sub-window)"
echo "[side-by-side] rust_full_interval=${RUST_FULL_INTERVAL} (Rust single active region)"
echo "[side-by-side] out=${OUT}"
echo "[side-by-side] compiling Java HcFullParityGateDump..."
"${JAVA_DUMP}" apply-summary "${REF}" "${BAM}" "${INT}" >/dev/null || true
./scripts/parity/run_hc_full_parity_java_compile.sh

run_pair() {
  local name="$1"
  shift
  local rust_cmd=("$@")
  echo "=== ${name} ==="
  "${RUST_DUMP}" "${rust_cmd[@]}" >"${OUT}/rust_${name}.tsv" 2>"${OUT}/rust_${name}.err" || true
  "${JAVA_DUMP}" "${rust_cmd[@]}" >"${OUT}/java_${name}.tsv" 2>"${OUT}/java_${name}.err" || true
  if diff -u "${OUT}/rust_${name}.tsv" "${OUT}/java_${name}.tsv" >"${OUT}/diff_${name}.txt" 2>&1; then
    echo "[${name}] MATCH"
  else
    echo "[${name}] DIFF (see ${OUT}/diff_${name}.txt)"
  fi
}

"${RUST_DUMP}" assembly-regions "${REF}" "${BAM}" "${RUST_FULL_INTERVAL}" >"${OUT}/rust_assembly-regions.tsv" 2>/dev/null || true
"${JAVA_DUMP}" assembly-regions "${REF}" "${BAM}" "${RUST_FULL_INTERVAL}" >"${OUT}/java_assembly-regions.tsv" 2>/dev/null || true
echo "=== assembly-regions (full interval ${RUST_FULL_INTERVAL}) ==="
cat "${OUT}/rust_assembly-regions.tsv" 2>/dev/null || true
cat "${OUT}/java_assembly-regions.tsv" 2>/dev/null || true

run_pair raw-activity raw-activity "${REF}" "${BAM}" "2:92307220-92307320" 100
run_pair apply-summary apply-summary "${REF}" "${BAM}" "${INT}"
run_pair assembly-region-haplotypes assembly-region-haplotypes "${REF}" "${BAM}" "${INT}"
run_pair assembly-region-kmer-probe assembly-region-kmer-probe "${REF}" "${BAM}" "${INT}"
run_pair assembly-region-assembly-stages assembly-region-assembly-stages "${REF}" "${BAM}" "${INT}"
run_pair assembly-region-assembly-stages-finalize assembly-region-assembly-stages-finalize "${REF}" "${BAM}" "${INT}"
run_pair assembly-region-finalize-reads assembly-region-finalize-reads "${REF}" "${BAM}" "${INT}"
run_pair assembly-region-kbest-paths assembly-region-kbest-paths "${REF}" "${BAM}" "${INT}"

echo ""
echo "--- haplotypes (first 15 data lines) ---"
echo "java:"
grep -v '^#' "${OUT}/java_assembly-region-haplotypes.tsv" 2>/dev/null | head -15 || true
echo "rust:"
grep -v '^#' "${OUT}/rust_assembly-region-haplotypes.tsv" 2>/dev/null | head -15 || true

echo ""
echo "--- k=85 probe row (expanded 85) ---"
awk -F'\t' '$2==85 && $1=="expanded"' "${OUT}/java_assembly-region-kmer-probe.tsv" 2>/dev/null | head -1 || true
awk -F'\t' '$2==85 && $1=="expanded"' "${OUT}/rust_assembly-region-kmer-probe.tsv" 2>/dev/null | head -1 || true

echo ""
echo "--- activity spike locus 92307296 (raw) ---"
awk -F'\t' 'NR==1 || $2==92307296' "${OUT}/rust_raw-activity.tsv" 2>/dev/null || \
  "${RUST_DUMP}" raw-activity "${REF}" "${BAM}" "2:92307296-92307296" 100 2>/dev/null | tee "${OUT}/rust_raw-activity-spike.tsv"
awk -F'\t' 'NR==1 || $2==92307296' "${OUT}/java_raw-activity.tsv" 2>/dev/null || \
  "${JAVA_DUMP}" raw-activity "${REF}" "${BAM}" "2:92307296-92307296" 100 2>/dev/null | tee "${OUT}/java_raw-activity-spike.tsv"

echo ""
echo "--- assembly stages (rt after_build / after_prune_dangling) ---"
awk -F'\t' '$1=="rt" && ($2 ~ /threading_after_build|threading_after_prune_dangling/)' "${OUT}/java_assembly-region-assembly-stages.tsv" 2>/dev/null || true
awk -F'\t' '$1=="rt" && ($2 ~ /threading_after_build|threading_after_prune_dangling/)' "${OUT}/rust_assembly-region-assembly-stages.tsv" 2>/dev/null || true

echo ""
echo "Full outputs: ${OUT}/"
