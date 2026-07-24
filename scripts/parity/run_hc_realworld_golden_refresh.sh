#!/usr/bin/env bash
# Regenerate committed L3 golden VCF (p11_java_positive / chrLive:1-63) and optional Java cross-check.
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

ref="${PARITY_HC_REALWORLD_REF:-${repo_root}/parity/fixtures/p5_live_reference.fa}"
sam="${repo_root}/parity/fixtures/p11_java_positive.sam"
cache="${repo_root}/parity/build/sam-indexed-bam/p11_java_positive.bam"
interval="${PARITY_HC_REALWORLD_INTERVAL:-chrLive:1-63}"
golden="${repo_root}/parity/fixtures/hc-full-parity/e2e-real/expected/p11_java_positive_chrlive_golden.vcf"
rust_out="${repo_root}/parity/reports/hc-realworld-golden-refresh.vcf"
java_out="${repo_root}/parity/reports/hc-realworld-golden-refresh-java.vcf"

mkdir -p "$(dirname "${golden}")" "$(dirname "${rust_out}")" "$(dirname "${cache}")"
if [[ ! -f "${cache}" ]]; then
  samtools view -bS "${sam}" | samtools sort -o "${cache}"
  samtools index "${cache}"
fi

echo "[hc-realworld-golden-refresh] Rust HC → ${rust_out}"
cargo run -q --bin gatk-rs -- \
  HaplotypeCaller \
  -R "${ref}" \
  -I "${cache}" \
  -L "${interval}" \
  -O "${rust_out}"

python3 "${repo_root}/scripts/parity/compare_hc_realworld_vcfs.py" \
  "${rust_out}" \
  --require-non-vacuous

# Normalize reference header to repo-relative path for committed golden.
sed 's|^##reference=.*|##reference=parity/fixtures/p5_live_reference.fa|' "${rust_out}" >"${golden}"

echo "[hc-realworld-golden-refresh] wrote ${golden}"

if [[ "${PARITY_SKIP_JAVA_REALWORLD_GOLDEN:-0}" != "1" ]]; then
  # shellcheck disable=SC1091
  source "${repo_root}/docs/GATK_PINNED.env"
  echo "[hc-realworld-golden-refresh] Java HC → ${java_out}"
  docker run --rm --platform "${GATK_DOCKER_PLATFORM}" \
    -v "${repo_root}:${repo_root}" \
    -w "${repo_root}" \
    "${GATK_DOCKER_IMAGE}" \
    gatk HaplotypeCaller \
    -R "${ref}" \
    -I "${cache}" \
    -L "${interval}" \
    -O "${java_out}" \
    --verbosity ERROR

  python3 "${repo_root}/scripts/parity/compare_hc_realworld_vcfs.py" \
    "${rust_out}" \
    --java "${java_out}" \
    --require-non-vacuous \
    --require-java-identity \
    --require-java-l3
fi

echo "[hc-realworld-golden-refresh] done"
