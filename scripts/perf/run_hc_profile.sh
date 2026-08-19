#!/usr/bin/env bash
# Production HC profile runner (observe-only).
#
# Usage:
#   ./scripts/perf/run_hc_profile.sh \
#     -R parity/realworld/assets/hs37d5.simple.fa \
#     -I parity/realworld/na12878_ci_loser_windows/01_chr21_w09.bam \
#     -L 21:9500000-9700000 \
#     -O /tmp/w09_profile.vcf \
#     --out-dir docs/perf/runs/hc_profile_w09
#
# Extra args after -- are passed to HaplotypeCaller.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

OUT_DIR="docs/perf/runs/hc_profile_latest"
THREADS="${RAYON_NUM_THREADS:-2}"
PAIR_HMM="${PAIR_HMM:-FASTEST_AVAILABLE}"
BIN="${GATK_RS_BIN:-target/release/gatk-rs}"
HC_ARGS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --bin) BIN="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --pair-hmm) PAIR_HMM="$2"; shift 2 ;;
    --) shift; HC_ARGS+=("$@"); break ;;
    *) HC_ARGS+=("$1"); shift ;;
  esac
done

if [[ ${#HC_ARGS[@]} -eq 0 ]]; then
  echo "error: pass HaplotypeCaller args (at least -R -I -L -O)" >&2
  exit 2
fi

mkdir -p "$OUT_DIR"
PROFILE_JSON="$OUT_DIR/hc_profile.json"
TRACE="$OUT_DIR/hc_rss.trace"

if [[ ! -x "$BIN" ]]; then
  echo "[profile] building release gatk-rs..."
  cargo build -p gatk-cli --release --locked
  BIN=target/release/gatk-rs
fi

export RAYON_NUM_THREADS="$THREADS"
# Pair with RSS TRACE for fine-grained assemble phase dual-write into the profiler.
export GATK_RS_HC_RSS_TRACE=1
export GATK_RS_HC_PROFILE="$PROFILE_JSON"
# Product wall shape (not Peak sequential).
unset GATK_RS_HC_SEQUENTIAL || true

echo "[profile] out=$OUT_DIR threads=$THREADS pair-hmm=$PAIR_HMM"
echo "[profile] profile=$PROFILE_JSON"
set +e
"$BIN" HaplotypeCaller \
  --threads "$THREADS" \
  --pair-hmm "$PAIR_HMM" \
  "${HC_ARGS[@]}" \
  >"$OUT_DIR/hc.stdout" 2>"$TRACE"
rc=$?
set -e

if [[ -f "$PROFILE_JSON" ]]; then
  echo "[profile] wrote $PROFILE_JSON"
  echo "[profile] wrote ${PROFILE_JSON%.json}.md"
else
  echo "[profile] WARNING: no profile JSON (did the process init HC_PROFILE?)" >&2
fi

if [[ -f "$TRACE" ]]; then
  python3 scripts/parity/giab/summarize_hc_rss_trace_wall.py "$TRACE" \
    | tee "$OUT_DIR/trace_summary.txt" || true
fi

exit "$rc"
