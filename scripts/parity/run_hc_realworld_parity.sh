#!/usr/bin/env bash
# L3 real-world HC parity — assembly-region-v1 vs golden VCF + Java CHROM/POS/REF/ALT (J.2.2 / CI.2).
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

golden_default="${repo_root}/parity/fixtures/hc-full-parity/e2e-real/expected/p11_java_positive_chrlive_golden.vcf"
if [[ -z "${PARITY_HC_REALWORLD_GOLDEN_VCF:-}" && -f "${golden_default}" ]]; then
  export PARITY_HC_REALWORLD_GOLDEN_VCF="${golden_default}"
fi

# Non-vacuous bundled corpus: p11_java_positive (alt haplotype, gq≥99) on p5_live reference.
if [[ -z "${PARITY_HC_REALWORLD_REF:-}" ]]; then
  export PARITY_HC_REALWORLD_REF="${repo_root}/parity/fixtures/p5_live_reference.fa"
fi
if [[ -z "${PARITY_HC_REALWORLD_BAM:-}" ]]; then
  sam="${repo_root}/parity/fixtures/p11_java_positive.sam"
  cache="${repo_root}/parity/build/sam-indexed-bam/p11_java_positive.bam"
  if [[ ! -f "${cache}" ]]; then
    mkdir -p "$(dirname "${cache}")"
    samtools view -bS "${sam}" | samtools sort -o "${cache}"
    samtools index "${cache}"
  fi
  export PARITY_HC_REALWORLD_BAM="${cache}"
fi
export PARITY_HC_REALWORLD_INTERVAL="${PARITY_HC_REALWORLD_INTERVAL:-chrLive:1-63}"

if [[ "${PARITY_HC_REALWORLD_STRICT:-0}" == "1" && -z "${PARITY_HC_REALWORLD_GOLDEN_VCF:-}" ]]; then
  echo "[hc-realworld-parity] PARITY_HC_REALWORLD_STRICT=1 requires PARITY_HC_REALWORLD_GOLDEN_VCF" >&2
  exit 1
fi

echo "[hc-realworld-parity] ref=${PARITY_HC_REALWORLD_REF}"
echo "[hc-realworld-parity] bam=${PARITY_HC_REALWORLD_BAM}"
echo "[hc-realworld-parity] interval=${PARITY_HC_REALWORLD_INTERVAL}"
echo "[hc-realworld-parity] golden=${PARITY_HC_REALWORLD_GOLDEN_VCF:-<none>} strict=${PARITY_HC_REALWORLD_STRICT:-0}"

"${repo_root}/scripts/parity/run_hc_full_parity_j_realworld.sh"

# L2 call-region VCF row (requires java_dumps from java_refresh).
if [[ "${PARITY_SKIP_HC_FULL_L2:-0}" != "1" ]]; then
  echo "[hc-realworld-parity] run_hc_full_parity_j2_vcf.sh (incl. call-region)"
  "${repo_root}/scripts/parity/run_hc_full_parity_j2_vcf.sh"
fi

echo "[hc-realworld-parity] done"
