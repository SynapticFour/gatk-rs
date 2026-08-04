#!/usr/bin/env bash
# Peak-RSS regression gate for HaplotypeCaller (PR Check).
#
# Catches:
#   A) SequenceDictionary::from_fasta_path materialising multi-Gb FASTA bodies
#   B) Unbounded k-best frontier on bushy/cyclic assembly graphs (NA12878
#      chr20:10098169 / 10098500 — the ~60 GiB realistic-window failure mode)
#   C) Contig-scale Smith-Waterman / PairHMM DP matrices
#
# Always (no genomic assets required):
#   1) Unit: dictionary prefers .fai
#   2) Unit: k-best heap/path bounds stay tight
#   3) Unit: SW refuses contig-scale matrices
#   4) dict_load_probe Peak-RSS on a synthetic fai-only FASTA
#
# When the NA12878 staging BAM + hs37d5.simple.fa are present:
#   5) HC on the minimal bomb window 20:10098500-10099500 (--threads 1);
#      Peak-RSS must stay under HC_RSS_MAX_MIB (default 256).
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"

# Avoid Cursor/sandbox alternate target dirs silently testing a stale binary.
export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-${repo_root}/target}"

HC_RSS_MAX_MIB="${HC_RSS_MAX_MIB:-256}"
DICT_RSS_MAX_MIB="${DICT_RSS_MAX_MIB:-64}"

echo "[hc-rss] unit: shared_bam Arc share + COW"
cargo test -p gatk-haplotypecaller --lib shared_bam -- --nocapture

echo "[hc-rss] unit: dictionary_prefers_fai_without_reading_fasta_body"
cargo test -p gatk-core --lib dictionary_prefers_fai_without_reading_fasta_body -- --nocapture

echo "[hc-rss] unit: kbest_bounds_are_finite_and_tight"
cargo test -p gatk-haplotypecaller --lib kbest_bounds_are_finite_and_tight -- --nocapture

echo "[hc-rss] unit: smith_waterman oversized matrix is refused"
cargo test -p gatk-haplotypecaller --lib oversized_matrix_is_refused -- --nocapture

echo "[hc-rss] unit: oversized assembly region read count is refused"
cargo test -p gatk-haplotypecaller --lib oversized_assembly_region_read_count_is_refused -- --nocapture

echo "[hc-rss] building dict_load_probe + gatk-rs (release)"
cargo build -p gatk-core --release --example dict_load_probe --locked
cargo build -p gatk-cli --release --locked --bin gatk-rs

tmp="$(mktemp -d)"
trap 'rm -rf "${tmp}"' EXIT

# Synthetic "genome-scale" index without a genome-scale FASTA body.
fa="${tmp}/poison.fa"
fai="${fa}.fai"
printf '>20\n' >"${fa}"
# name length offset linebases linewidth — LN matches hs37d5 chr20; body intentionally empty.
printf '20\t63025520\t4\t60\t61\n' >"${fai}"

parse_rss_mib() {
  local time_log="$1"
  python3 - "$time_log" <<'PY'
import re, sys
text = open(sys.argv[1], encoding="utf-8", errors="replace").read()
# GNU time -v
m = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", text)
if m:
    print(f"{int(m.group(1)) / 1024.0:.3f}")
    raise SystemExit(0)
# macOS /usr/bin/time -l (bytes before label; locale decimal may use comma elsewhere)
m = re.search(r"(\d+)\s+maximum resident set size\b", text)
if m:
    print(f"{int(m.group(1)) / (1024.0 * 1024.0):.3f}")
    raise SystemExit(0)
print("nan")
raise SystemExit(2)
PY
}

run_timed() {
  local log="$1"
  shift
  if /usr/bin/time -v true >/dev/null 2>&1; then
    /usr/bin/time -v -o "${log}" "$@"
  elif command -v gtime >/dev/null 2>&1 && gtime -v true >/dev/null 2>&1; then
    gtime -v -o "${log}" "$@"
  else
    /usr/bin/time -l "$@" >"${log}.stdout" 2>"${log}"
  fi
}

echo "[hc-rss] dict_load_probe on fai-only FASTA (max ${DICT_RSS_MAX_MIB} MiB)"
dict_time="${tmp}/dict.time"
run_timed "${dict_time}" \
  "${CARGO_TARGET_DIR}/release/examples/dict_load_probe" "${fa}"
dict_mib="$(parse_rss_mib "${dict_time}")"
echo "[hc-rss] dict Peak-RSS=${dict_mib} MiB"
python3 -c "import sys; m=float('${dict_mib}'); sys.exit(0 if m <= float('${DICT_RSS_MAX_MIB}') else 1)" \
  || {
    echo "[hc-rss] FAIL: dict Peak-RSS ${dict_mib} MiB > ${DICT_RSS_MAX_MIB} MiB" >&2
    exit 1
  }

REAL_REF="${HC_RSS_REF:-${repo_root}/parity/realworld/assets/hs37d5.simple.fa}"
REAL_BAM="${HC_RSS_BAM:-${repo_root}/parity/realworld/na12878_giab_window_mem_2mb_b37/NA12878_giab_window.b37.bam}"
# Minimal window that previously OOM'd assemble/HC via unbounded k-best (not the 50 kb gate).
REAL_INTERVAL="${HC_RSS_INTERVAL:-20:10098500-10099500}"

# Accept either `foo.bam.bai` or Picard-style `foo.bai` beside `foo.bam`.
real_bai=""
if [[ -f "${REAL_BAM}.bai" ]]; then
  real_bai="${REAL_BAM}.bai"
elif [[ -f "${REAL_BAM%.bam}.bai" ]]; then
  real_bai="${REAL_BAM%.bam}.bai"
fi

if [[ -f "${REAL_REF}" && -f "${REAL_REF}.fai" && -f "${REAL_BAM}" && -n "${real_bai}" ]]; then
  echo "[hc-rss] HC ${REAL_INTERVAL} threads=1 (max ${HC_RSS_MAX_MIB} MiB)"
  hc_time="${tmp}/hc.time"
  hc_vcf="${tmp}/hc.vcf"
  set +e
  run_timed "${hc_time}" \
    "${CARGO_TARGET_DIR}/release/gatk-rs" HaplotypeCaller \
      -R "${REAL_REF}" -I "${REAL_BAM}" -O "${hc_vcf}" \
      -L "${REAL_INTERVAL}" --threads 1
  hc_rc=$?
  set -e
  if [[ "${hc_rc}" -ne 0 ]]; then
    echo "[hc-rss] FAIL: HaplotypeCaller exited ${hc_rc}" >&2
    tail -40 "${hc_time}" >&2 || true
    exit 1
  fi
  hc_mib="$(parse_rss_mib "${hc_time}")"
  echo "[hc-rss] HC Peak-RSS=${hc_mib} MiB"
  python3 -c "import sys; m=float('${hc_mib}'); sys.exit(0 if m <= float('${HC_RSS_MAX_MIB}') else 1)" \
    || {
      echo "[hc-rss] FAIL: HC Peak-RSS ${hc_mib} MiB > ${HC_RSS_MAX_MIB} MiB" >&2
      exit 1
    }
else
  echo "[hc-rss] skip HC window (REALISTIC BAM/ref not staged on this runner)"
fi

echo "[hc-rss] OK"
