#!/usr/bin/env bash
# Fail unless pinned Java GATK 4.4 loads a native vector PairHMM for FASTEST_AVAILABLE.
# See docs/ci/PERF_BENCHMARK_HOST.md and docs/GATK_PINNED.env
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"
# shellcheck source=../parity/lib_pinned_gatk.sh
source "${repo_root}/scripts/parity/lib_pinned_gatk.sh"

REF="${HC_MEM_REF:-${repo_root}/parity/fixtures/reference.fa}"
BAM="${HC_MEM_BAM:-${repo_root}/parity/fixtures/sample.bam}"
INTERVAL="${HC_MEM_INTERVAL:-chr1:1-32}"
outdir="$(mktemp -d "${TMPDIR:-/tmp}/gatk-avx-verify.XXXXXX")"
trap 'rm -rf "${outdir}"' EXIT
out_vcf="${outdir}/verify.vcf"
log="${outdir}/gatk.log"

# Prefer ripgrep; fall back to grep -E (self-hosted images may lack rg).
grep_i() {
  local pat="$1" file="$2"
  if command -v rg >/dev/null 2>&1; then
    rg -qi "${pat}" "${file}"
  else
    grep -Eiq "${pat}" "${file}"
  fi
}

grep_ni() {
  local pat="$1" file="$2"
  if command -v rg >/dev/null 2>&1; then
    rg -ni "${pat}" "${file}"
  else
    grep -Ein "${pat}" "${file}" || true
  fi
}

persist_log() {
  if [[ -n "${PERF_RUN_DIR:-}" && -f "${log}" ]]; then
    mkdir -p "${PERF_RUN_DIR}"
    cp -f "${log}" "${PERF_RUN_DIR}/java_avx_verify.log"
  fi
}

if [[ ! -f "${BAM}.bai" ]]; then
  samtools index "${BAM}"
fi

echo "[avx-verify] image=${GATK_DOCKER_IMAGE} sha=${GATK_PINNED_SHA}"
echo "[avx-verify] interval=${INTERVAL}"

host_root="${repo_root}"
plat_args=()
if [[ "$(uname -s)" == "Darwin" ]]; then
  plat_args+=(--platform "${GATK_DOCKER_PLATFORM:-linux/amd64}")
fi

docker run --rm \
  "${plat_args[@]}" \
  -v "${host_root}:${host_root}" \
  -w "${host_root}" \
  "${GATK_DOCKER_IMAGE}" \
  gatk --java-options '-Xms1g -Xmx2g' HaplotypeCaller \
    -R "${REF}" -I "${BAM}" -O "${out_vcf}" -L "${INTERVAL}" \
    --pair-hmm-implementation FASTEST_AVAILABLE \
  >"${log}" 2>&1 || {
    persist_log
    echo "[avx-verify] GATK failed; last 80 log lines:" >&2
    tail -80 "${log}" >&2 || true
    exit 1
  }

persist_log

# Broad GATK / GKL typically logs the selected implementation. Accept common phrases.
if grep_i 'Using PairHMM implementation:[[:space:]]*(AVX|OMP|LOGLESS_AVX)|VectorLoglessPairHMM|NativePairHMM|IntelPairHmm|AVXPairHMM' "${log}"; then
  echo "[avx-verify] OK — native / vector PairHMM path detected:"
  grep_ni 'PairHMM|AVX|VectorLogless|IntelPairHmm|NativePairHMM' "${log}" | head -20
  exit 0
fi

# Some images only log at DEBUG; also accept GKL load success lines.
if grep_i 'libgkl_pairhmm|GKL.*PairHMM.*(loaded|initialized)|Successfully loaded.*pairhmm' "${log}"; then
  echo "[avx-verify] OK — GKL PairHMM library load detected:"
  grep_ni 'gkl|pairhmm|avx' "${log}" | head -20
  exit 0
fi

echo "[avx-verify] FAIL — could not confirm native AVX PairHMM; possible Java fallback." >&2
echo "[avx-verify] Relevant log lines:" >&2
grep_ni 'PairHMM|AVX|LOGLESS|LOG10|implementation|GKL|libgkl' "${log}" | head -40 >&2 || true
echo "[avx-verify] Full log saved under TEMP (and PERF_RUN_DIR if set)" >&2
exit 1
