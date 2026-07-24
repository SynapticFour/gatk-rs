#!/usr/bin/env bash
# Diagnose P12 cluster assembly (92307200–92307363) without cargo — uses prebuilt release dump binary.
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"

REF="${P12_REFERENCE:-${repo_root}/parity/realworld/assets/hs37d5.simple.fa}"
BAM="${P12_BAM:-${repo_root}/parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam}"
INT="${P12_CLUSTER_INTERVAL:-2:92307200-92307400}"
DUMP="${PARITY_DUMP_BIN:-${repo_root}/target/release/examples/hc_full_parity_gate_dump}"

if [[ ! -x "${DUMP}" ]]; then
  echo "Build dump example first: ./scripts/parity/build_gatk_rs_release.sh && \\"
  echo "  CARGO_TARGET_DIR=${repo_root}/target cargo build --release -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump"
  exit 1
fi

echo "[diagnose] interval=${INT}"
echo "[diagnose] Java side-by-side: ./scripts/parity/diagnose_p12_cluster_assembly_java_rust.sh"
echo "=== apply-summary ==="
"${DUMP}" apply-summary "${REF}" "${BAM}" "${INT}"
echo "=== assembly-region-haplotypes (first active) ==="
"${DUMP}" assembly-region-haplotypes "${REF}" "${BAM}" "${INT}" active | head -20
echo "=== assembly-region-kmer-probe (first active) ==="
"${DUMP}" assembly-region-kmer-probe "${REF}" "${BAM}" "${INT}"
echo "=== assembly-region-assembly-stages k=85 (first active) ==="
"${DUMP}" assembly-region-assembly-stages "${REF}" "${BAM}" "${INT}"
echo "=== reference slice 92307320-92307335 ==="
samtools faidx "${REF}" 2:92307320-92307335 2>/dev/null || true
rust_vcf="${repo_root}/parity/reports/p12_realworld_na12878_20k.rust.vcf"
java_vcf="${repo_root}/parity/reports/p12_realworld_na12878_20k.java.vcf"
if [[ -s "${rust_vcf}" && -s "${java_vcf}" ]]; then
  echo "=== cluster VCF rows (92307320-92307340) ==="
  awk -v s=92307320 -v e=92307340 '$1=="2" && $2>=s && $2<=e {print "rust",$0}' "${rust_vcf}" | grep -v '^#' || true
  awk -v s=92307320 -v e=92307340 '$1=="2" && $2>=s && $2<=e {print "java",$0}' "${java_vcf}" | grep -v '^#' || true
  json_tmp="$(mktemp)"
  python3 "${repo_root}/scripts/parity/p12_na12878_summarize.py" \
    "${json_tmp}" "${json_tmp}.md" "${java_vcf}" "${rust_vcf}" 0 0 '{}'
  rm -f "${json_tmp}" "${json_tmp}.md"
fi
