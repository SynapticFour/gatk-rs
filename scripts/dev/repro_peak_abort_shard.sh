#!/usr/bin/env bash
# Reproducible Peak / RSS-abort probe for a GIAB-style HC shard.
#
# ## Proven hotspot (2026-08-11, tip with abort watchdog)
#
# CI shard `00_chr20_w20` (`20:20000001-21000000`) hard-OOMs because Peak climbs
# inside `build_threading_graph_core` — *before* bushiness skip / k-best abort:
#
#   locus  20:20301092-20301239  (reads=112)
#   phase  rt_graph_build_begin kmer=35 before_remove_paths=true
#   signal HC_RSS_ABORT_WATCHDOG fires; HC_RSS_ABORT (k-best soft-land) does NOT
#
# Minimal local repro (BAM must cover the locus; w20 slice works):
#   SKIP_SLICE=1 INTERVAL=20:20300000-20302000 ABORT_MIB=256 \
#     BAM=parity/giab/runs/local-peak-repro/HG001.20-20000001-21000000.bam \
#     ./scripts/dev/repro_peak_abort_shard.sh
#
# Full CI window (slow; same failure class):
#   INTERVAL=20:20000001-21000000 ABORT_MIB=512 ./scripts/dev/repro_peak_abort_shard.sh
#
# Named historical spike (now soft; ~40 MiB — not the CI failure):
#   INTERVAL=20:10098000-10098600 SKIP_SLICE=1 \
#   BAM=parity/realworld/na12878_giab_window_mem_500kb_b37/NA12878_giab_window.b37.bam \
#   ./scripts/dev/repro_peak_abort_shard.sh
#
# Log signals:
#   HC_RSS_ABORT_CONFIG   — in-process limit was parsed
#   HC_RSS_TRACE phase=…  — which phase/locus was active as RSS climbed
#   HC_RSS_ABORT_WATCHDOG — sampler saw RSS ≥ limit (even if k-best never checked)
#   HC_RSS_ABORT          — soft-abort path actually ran (k-best / ingest check)
#
# Usage:
#   ./scripts/dev/repro_peak_abort_shard.sh
#   ABORT_MIB=512 INTERVAL=20:20000001-21000000 ./scripts/dev/repro_peak_abort_shard.sh
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-${repo_root}/target}"

# shellcheck source=../parity/giab/lib_giab.sh
source "${repo_root}/scripts/parity/giab/lib_giab.sh"

INTERVAL="${INTERVAL:-20:20000001-21000000}"
ABORT_MIB="${ABORT_MIB:-512}"
THREADS="${THREADS:-2}"
OUT_DIR="${OUT_DIR:-${repo_root}/parity/giab/runs/local-peak-repro}"
REF="${REF:-${repo_root}/parity/realworld/assets/hs37d5.simple.fa}"
SKIP_SLICE="${SKIP_SLICE:-0}"
ftp_data="${GIAB_FTP_DATA:-https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data}"
BAM_URL="${GIAB_HG001_BAM_URL:-${ftp_data}/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam}"

mkdir -p "${OUT_DIR}"
log="${OUT_DIR}/hc-${INTERVAL//[:]/-}.log"
vcf="${OUT_DIR}/hc-${INTERVAL//[:]/-}.vcf"
bam="${BAM:-${OUT_DIR}/HG001.${INTERVAL//[:]/-}.bam}"

if [[ ! -f "${REF}" ]]; then
  echo "missing REF=${REF}" >&2
  exit 1
fi

if [[ "${SKIP_SLICE}" != "1" ]]; then
  if [[ ! -f "${bam}" || ! -f "${bam}.bai" ]]; then
    echo "[repro] slicing ${INTERVAL} from ${BAM_URL}" >&2
    bai_local="${OUT_DIR}/remote.bai"
    if [[ ! -f "${bai_local}" ]]; then
      curl -fL --retry 3 -o "${bai_local}.partial" "${BAM_URL}.bai"
      mv -f "${bai_local}.partial" "${bai_local}"
    fi
    samtools view -b -X "${BAM_URL}" "${bai_local}" "${INTERVAL}" > "${bam}.partial"
    mv -f "${bam}.partial" "${bam}"
    samtools index "${bam}" "${bam}.bai"
  else
    echo "[repro] reuse BAM ${bam}" >&2
  fi
else
  if [[ ! -f "${bam}" ]]; then
    echo "SKIP_SLICE=1 but BAM missing: ${bam}" >&2
    exit 1
  fi
fi

echo "[repro] building release gatk-rs" >&2
cargo build -p gatk-cli --release --locked --bin gatk-rs

bin="${CARGO_TARGET_DIR}/release/gatk-rs"
echo "[repro] ABORT_MIB=${ABORT_MIB} INTERVAL=${INTERVAL} THREADS=${THREADS}" >&2
echo "[repro] log → ${log}" >&2

set +e
env \
  GATK_RS_HC_SEQUENTIAL=1 \
  GATK_RS_HC_RSS_ABORT_MIB="${ABORT_MIB}" \
  GATK_RS_HC_RSS_TRACE=1 \
  RAYON_NUM_THREADS="${THREADS}" \
  MALLOC_ARENA_MAX=2 \
  "${bin}" HaplotypeCaller \
    -R "${REF}" \
    -I "${bam}" \
    -O "${vcf}" \
    -L "${INTERVAL}" \
  >"${log}.stdout" 2>"${log}"
rc=$?
set -e

echo "[repro] exit=${rc}" >&2
echo "[repro] --- abort / phase signals ---" >&2
rg -n "HC_RSS_ABORT|HC_RSS_TRACE phase=|HC_RSS_TRACE sample" "${log}" | head -80 || true
echo "[repro] --- peak sample ---" >&2
rg -n "HC_RSS_TRACE sample" "${log}" | tail -5 || true
exit "${rc}"
