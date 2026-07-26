#!/usr/bin/env bash
# J-D01 / J-D03 / J.2.2 — real-world VCF parity (p11 non-vacuous corpus + golden VCF).
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"

if [[ -z "${PARITY_HC_REALWORLD_BAM:-}" || -z "${PARITY_HC_REALWORLD_REF:-}" ]]; then
  echo "[hc-full-parity-j-realworld] skipped (set PARITY_HC_REALWORLD_BAM and PARITY_HC_REALWORLD_REF)"
  exit 0
fi

export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"
interval="${PARITY_HC_REALWORLD_INTERVAL:-chrLive:1-63}"
out_vcf="${repo_root}/parity/reports/hc-realworld-out.vcf"
golden_default="${repo_root}/parity/fixtures/hc-full-parity/e2e-real/expected/p11_java_positive_chrlive_golden.vcf"
golden="${PARITY_HC_REALWORLD_GOLDEN_VCF:-${golden_default}}"
java_out="${repo_root}/parity/reports/hc-realworld-java.vcf"
strict="${PARITY_HC_REALWORLD_STRICT:-0}"

if [[ "${strict}" == "1" && ! -f "${golden}" ]]; then
  echo "[hc-full-parity-j-realworld] PARITY_HC_REALWORLD_STRICT=1 requires golden at ${golden}" >&2
  exit 1
fi

mkdir -p "${repo_root}/parity/reports"

echo "[hc-full-parity-j-realworld] running gatk-rs HC on ${PARITY_HC_REALWORLD_BAM} interval=${interval}"
cargo run -q --bin gatk-rs -- \
  HaplotypeCaller \
  -R "${PARITY_HC_REALWORLD_REF}" \
  -I "${PARITY_HC_REALWORLD_BAM}" \
  -L "${interval}" \
  -O "${out_vcf}" \
  2>"${repo_root}/parity/reports/hc-realworld.stderr"

if ! grep -q 'assembly-region-v1' "${out_vcf}"; then
  echo "[hc-full-parity-j-realworld] expected GATK_RS_HC_PIPELINE=assembly-region-v1" >&2
  grep 'GATK_RS_HC_PIPELINE' "${out_vcf}" >&2 || true
  exit 1
fi

compare_args=(
  "${repo_root}/scripts/parity/compare_hc_realworld_vcfs.py"
  "${out_vcf}"
  --require-non-vacuous
)

if [[ -f "${golden}" ]]; then
  compare_args+=(--golden "${golden}")
fi

if [[ "${PARITY_HC_REALWORLD_JAVA_CHECK:-1}" == "1" && "${PARITY_SKIP_JAVA_REALWORLD_GOLDEN:-0}" != "1" ]]; then
  # shellcheck disable=SC1091
  source "${repo_root}/docs/GATK_PINNED.env"
  echo "[hc-full-parity-j-realworld] Java HC identity check"
  docker run --rm --platform "${GATK_DOCKER_PLATFORM}" \
    -v "${repo_root}:${repo_root}" \
    -w "${repo_root}" \
    "${GATK_DOCKER_IMAGE}" \
    gatk HaplotypeCaller \
    -R "${PARITY_HC_REALWORLD_REF}" \
    -I "${PARITY_HC_REALWORLD_BAM}" \
    -L "${interval}" \
    -O "${java_out}" \
    --verbosity ERROR \
    2>"${repo_root}/parity/reports/hc-realworld-java.stderr"
  compare_args+=(--java "${java_out}")
  if [[ "${strict}" == "1" ]]; then
    compare_args+=(--require-java-identity --require-java-l3 --require-java-l4)
  fi
fi

python3 "${compare_args[@]}"

echo "[hc-full-parity-j-realworld] L3 real-world VCF parity OK"
