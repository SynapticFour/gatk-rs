#!/usr/bin/env bash
# ASM-1: Java finalizeRegion vs Rust production finalize (reads + assembly stages).
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"

REF="${P12_REFERENCE:-${repo_root}/parity/realworld/assets/hs37d5.simple.fa}"
BAM="${P12_BAM:-${repo_root}/parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam}"
INT="${P12_CLUSTER_INTERVAL:-2:92307228-92307400}"
OUT="${ASM_FINALIZE_OUT:-${repo_root}/parity/reports/asm_finalize_parity}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"
RUST_DUMP="${PARITY_DUMP_BIN:-${CARGO_TARGET_DIR}/release/examples/hc_full_parity_gate_dump}"
JAVA_DUMP="${repo_root}/scripts/parity/run_hc_full_parity_java_dump.sh"

mkdir -p "${OUT}"

if [[ ! -x "${RUST_DUMP}" ]]; then
  echo "Build Rust dump: CARGO_TARGET_DIR=${repo_root}/target cargo build --release -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump"
  exit 1
fi

echo "[asm-finalize] interval=${INT}"
echo "[asm-finalize] out=${OUT}"
"${JAVA_DUMP}" apply-summary "${REF}" "${BAM}" "${INT}" >/dev/null || true
./scripts/parity/run_hc_full_parity_java_compile.sh

run_pair() {
  local name="$1"
  shift
  echo "=== ${name} ==="
  "${RUST_DUMP}" "$@" >"${OUT}/rust_${name}.tsv" 2>"${OUT}/rust_${name}.err" || true
  "${JAVA_DUMP}" "$@" >"${OUT}/java_${name}.tsv" 2>"${OUT}/java_${name}.err" || true
  if diff -u "${OUT}/rust_${name}.tsv" "${OUT}/java_${name}.tsv" >"${OUT}/diff_${name}.txt" 2>&1; then
    echo "[${name}] MATCH"
  else
    echo "[${name}] DIFF (see ${OUT}/diff_${name}.txt)"
  fi
}

run_pair assembly-region-finalize-reads assembly-region-finalize-reads "${REF}" "${BAM}" "${INT}"
run_pair assembly-region-assembly-stages-finalize assembly-region-assembly-stages-finalize "${REF}" "${BAM}" "${INT}"

echo ""
echo "--- finalize reads (phase=finalize rows) ---"
awk -F'\t' 'NR==1 || ($1=="read" && $3=="finalize")' "${OUT}/java_assembly-region-finalize-reads.tsv" 2>/dev/null | head -20 || true
awk -F'\t' 'NR==1 || ($1=="read" && $3=="finalize")' "${OUT}/rust_assembly-region-finalize-reads.tsv" 2>/dev/null | head -20 || true

echo ""
echo "--- assembly stages finalize (rt prune_before / dangling) ---"
awk -F'\t' '$1=="rt" && ($2 ~ /threading_after_prune_before_dangling|threading_after_dangling/)' \
  "${OUT}/java_assembly-region-assembly-stages-finalize.tsv" 2>/dev/null || true
awk -F'\t' '$1=="rt" && ($2 ~ /threading_after_prune_before_dangling|threading_after_dangling/)' \
  "${OUT}/rust_assembly-region-assembly-stages-finalize.tsv" 2>/dev/null || true

echo ""
echo "--- materialize vs finalize (Java only) ---"
"${JAVA_DUMP}" assembly-region-assembly-stages "${REF}" "${BAM}" "${INT}" \
  >"${OUT}/java_assembly-region-assembly-stages-materialize.tsv" 2>/dev/null || true
awk -F'\t' '$1=="dangling_recovery" || ($1=="rt" && $2 ~ /threading_after_prune_before/)' \
  "${OUT}/java_assembly-region-assembly-stages-materialize.tsv" 2>/dev/null || true
awk -F'\t' '$1=="dangling_recovery" || ($1=="rt" && $2 ~ /threading_after_prune_before/)' \
  "${OUT}/java_assembly-region-assembly-stages-finalize.tsv" 2>/dev/null || true

echo ""
echo "Full outputs: ${OUT}/"
