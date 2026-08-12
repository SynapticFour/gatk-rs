#!/usr/bin/env bash
# Differential wall-time / Peak-RSS: Java GATK 4.4 vs gatk-rs on a pinned GIAB window.
#
# Motivation (ci-subset 31523602882): cancelled rust shards timed out at 360 min while
# Java finished the same 1 Mb windows in 3–10 min. TRACE attributed ~hours to the
# post-assemble PairHMM gap (default rust = LOG10_PAIRHMM scalar).
#
# Default pin = neighborhood of the worst w29 region (20:29455745-29455993):
#   INTERVAL=20:29450000-29460000
#
# Usage:
#   ./scripts/dev/repro_pairhmm_java_rust_diff.sh
#   INTERVAL=20:29455700-29456000 PAIR_HMM=FASTEST_AVAILABLE ./scripts/dev/repro_pairhmm_java_rust_diff.sh
#   SKIP_JAVA=1 PAIR_HMM=LOG10_PAIRHMM ./scripts/dev/repro_pairhmm_java_rust_diff.sh
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-${repo_root}/target}"

# shellcheck source=../parity/giab/lib_giab.sh
source "${repo_root}/scripts/parity/giab/lib_giab.sh"
# shellcheck source=../parity/lib_pinned_gatk.sh
source "${repo_root}/scripts/parity/lib_pinned_gatk.sh"

INTERVAL="${INTERVAL:-20:29450000-29460000}"
THREADS="${THREADS:-2}"
PAIR_HMM="${PAIR_HMM:-}" # empty = rust production default (LOG10_PAIRHMM)
SKIP_JAVA="${SKIP_JAVA:-0}"
SKIP_RUST="${SKIP_RUST:-0}"
SKIP_SLICE="${SKIP_SLICE:-0}"
OUT_DIR="${OUT_DIR:-${repo_root}/parity/giab/runs/local-pairhmm-diff}"
REF="${REF:-${repo_root}/parity/realworld/assets/hs37d5.simple.fa}"
ftp_data="${GIAB_FTP_DATA:-https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data}"
BAM_URL="${GIAB_HG001_BAM_URL:-${ftp_data}/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam}"

mkdir -p "${OUT_DIR}"
tag="${INTERVAL//[:]/-}"
bam="${BAM:-${OUT_DIR}/HG001.${tag}.bam}"
rust_log="${OUT_DIR}/rust.${tag}.log"
java_log="${OUT_DIR}/java.${tag}.log"
rust_vcf="${OUT_DIR}/rust.${tag}.vcf"
java_vcf="${OUT_DIR}/java.${tag}.vcf"
rust_time="${OUT_DIR}/rust.${tag}.time.txt"
java_time="${OUT_DIR}/java.${tag}.time.txt"

if [[ ! -f "${REF}" ]]; then
  echo "missing REF=${REF}" >&2
  exit 1
fi

if [[ "${SKIP_SLICE}" != "1" ]]; then
  if [[ ! -f "${bam}" || ! -f "${bam}.bai" ]]; then
    echo "[diff] slicing ${INTERVAL} from ${BAM_URL}" >&2
    bai_local="${OUT_DIR}/remote.bai"
    if [[ ! -f "${bai_local}" ]]; then
      curl -fL --retry 3 -o "${bai_local}.partial" "${BAM_URL}.bai"
      mv -f "${bai_local}.partial" "${bai_local}"
    fi
    samtools view -b -X "${BAM_URL}" "${bai_local}" "${INTERVAL}" >"${bam}.partial"
    mv -f "${bam}.partial" "${bam}"
    samtools index "${bam}" "${bam}.bai"
  else
    echo "[diff] reuse BAM ${bam}" >&2
  fi
else
  [[ -f "${bam}" ]] || {
    echo "SKIP_SLICE=1 but BAM missing: ${bam}" >&2
    exit 1
  }
fi

run_timed() {
  local time_log="$1"
  local out_log="$2"
  shift 2
  # Keep HC stderr (TRACE) in out_log; time accounting in time_log.
  if command -v gtime >/dev/null 2>&1; then
    gtime -v -o "${time_log}" "$@" >"${out_log}.stdout" 2>"${out_log}"
  else
    # macOS /usr/bin/time -l writes stats on stderr after the command's stderr.
    {
      /usr/bin/time -l "$@" 2> >(tee "${out_log}" >&2)
    } >"${out_log}.stdout" 2>"${time_log}"
  fi
}

parse_rss() {
  python3 - "$1" <<'PY'
import re, sys
text = open(sys.argv[1], encoding="utf-8", errors="replace").read()
m = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", text)
if m:
    print(f"{int(m.group(1))/1024:.1f}")
    raise SystemExit(0)
m = re.search(r"(\d+)\s+maximum resident set size\b", text)
if m:
    print(f"{int(m.group(1))/(1024*1024):.1f}")
    raise SystemExit(0)
print("?")
PY
}

summarize_rust_trace() {
  python3 - "$1" <<'PY'
import re, sys
from statistics import median
text = open(sys.argv[1], encoding="utf-8", errors="replace").read()
rows = []
for ln in text.splitlines():
    if "phase=after_pairhmm" not in ln:
        continue
    ms = re.search(r"pairhmm_ms=(\d+)", ln)
    haps = re.search(r"haps=(\d+)", ln)
    reads = re.search(r"geno_reads=(\d+)", ln)
    impl = re.search(r"impl=(\S+)", ln)
    backend = re.search(r"backend=(\S+)", ln)
    loc = re.search(r"locus=(\S+)", ln)
    if not ms:
        continue
    rows.append(
        (
            int(ms.group(1)),
            int(haps.group(1)) if haps else 0,
            int(reads.group(1)) if reads else 0,
            impl.group(1) if impl else "?",
            backend.group(1) if backend else "?",
            loc.group(1) if loc else "?",
        )
    )
grows = []
preps = []
assigns = []
for ln in text.splitlines():
    if "phase=after_genotype" in ln:
        ms = re.search(r"genotype_ms=(\d+)", ln)
        asn = re.search(r"assign_genotype_ms=(\d+)", ln)
        if ms:
            grows.append(int(ms.group(1)))
        if asn:
            assigns.append(int(asn.group(1)))
    if "phase=before_assign_genotype" in ln:
        ms = re.search(r"post_pairhmm_prep_ms=(\d+)", ln)
        if ms:
            preps.append(int(ms.group(1)))
if not rows:
    print("  (no after_pairhmm TRACE — rebuild gatk-rs with instrumentation)")
else:
    vals = [r[0] for r in rows]
    print(f"  regions_with_pairhmm={len(rows)}  impl={rows[0][3]} backend={rows[0][4]}")
    print(
        f"  pairhmm_ms: sum={sum(vals)/1000:.1f}s  mean={sum(vals)/len(vals):.0f}  "
        f"median={median(vals):.0f}  max={max(vals)}"
    )
    if preps:
        print(
            f"  post_pairhmm_prep_ms: sum={sum(preps)/1000:.1f}s  mean={sum(preps)/len(preps):.0f}  "
            f"max={max(preps)}"
        )
    if assigns:
        print(
            f"  assign_genotype_ms: sum={sum(assigns)/1000:.1f}s  mean={sum(assigns)/len(assigns):.0f}  "
            f"max={max(assigns)}"
        )
    if grows:
        print(
            f"  genotype_ms(total_post_pairhmm): sum={sum(grows)/1000:.1f}s  mean={sum(grows)/len(grows):.0f}  "
            f"max={max(grows)}"
        )
    rows.sort(reverse=True)
    print("  top pairhmm_ms:")
    for ms, haps, reads, impl, backend, loc in rows[:8]:
        print(f"    {ms:7}ms  haps={haps:3} reads={reads:4}  {loc}")
PY
}

if [[ "${SKIP_RUST}" != "1" ]]; then
  echo "[diff] building release gatk-rs" >&2
  cargo build -p gatk-cli --release --locked --bin gatk-rs
  bin="${CARGO_TARGET_DIR}/release/gatk-rs"
  rust_args=(HaplotypeCaller -R "${REF}" -I "${bam}" -O "${rust_vcf}" -L "${INTERVAL}")
  if [[ -n "${PAIR_HMM}" ]]; then
    rust_args+=(--pair-hmm "${PAIR_HMM}")
  fi
  echo "[diff] RUST INTERVAL=${INTERVAL} THREADS=${THREADS} PAIR_HMM=${PAIR_HMM:-default}" >&2
  set +e
  run_timed "${rust_time}" "${rust_log}" env \
    GATK_RS_HC_SEQUENTIAL=1 \
    GATK_RS_HC_RSS_ABORT_MIB="${ABORT_MIB:-8192}" \
    GATK_RS_HC_RSS_TRACE=1 \
    RAYON_NUM_THREADS="${THREADS}" \
    MALLOC_ARENA_MAX=2 \
    "${bin}" "${rust_args[@]}"
  rust_rc=$?
  set -e
  echo "[diff] rust exit=${rust_rc} PeakRSS_MiB=$(parse_rss "${rust_time}")" >&2
  echo "[diff] rust PairHMM TRACE summary:" >&2
  summarize_rust_trace "${rust_log}"
  # Fallback: TRACE may have landed in time log on older script versions.
  if ! rg -q "phase=after_pairhmm" "${rust_log}" 2>/dev/null; then
    summarize_rust_trace "${rust_time}"
  fi
fi

if [[ "${SKIP_JAVA}" != "1" ]]; then
  img="${GATK_DOCKER_IMAGE:-us.gcr.io/broad-gatk/gatk:4.4.0.0}"
  plat="${GATK_DOCKER_PLATFORM:-linux/amd64}"
  echo "[diff] JAVA docker ${img} INTERVAL=${INTERVAL} threads=${THREADS}" >&2
  set +e
  run_timed "${java_time}" "${java_log}" docker run --rm --platform "${plat}" \
    -v "${repo_root}:${repo_root}" -w "${repo_root}" \
    "${img}" gatk HaplotypeCaller \
    -R "${REF}" -I "${bam}" -O "${java_vcf}" -L "${INTERVAL}" \
    --verbosity ERROR \
    --native-pair-hmm-threads "${THREADS}"
  java_rc=$?
  set -e
  echo "[diff] java exit=${java_rc} PeakRSS_MiB=$(parse_rss "${java_time}")" >&2
  if [[ -f "${java_time}" ]]; then
    rg -n "Elapsed|real|User time|System time|maximum resident" "${java_time}" | head -20 || true
  fi
fi

echo "[diff] artifacts under ${OUT_DIR}" >&2
