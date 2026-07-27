#!/usr/bin/env bash
# Shared helpers for GIAB genome-wide equivalence runners.
# shellcheck shell=bash

giab_time_backend() {
  # Prefer GNU time -v (Linux CI / gtime on macOS); else macOS /usr/bin/time -l.
  if /usr/bin/time -v true >/dev/null 2>&1; then
    echo "gnu"
  elif command -v gtime >/dev/null 2>&1 && gtime -v true >/dev/null 2>&1; then
    echo "gtime"
  else
    echo "macos"
  fi
}

giab_run_timed() {
  # Args: backend logfile command...
  local backend="$1" log="$2"
  shift 2
  case "${backend}" in
    gnu) /usr/bin/time -v -o "${log}" "$@" ;;
    gtime) gtime -v -o "${log}" "$@" ;;
    macos) /usr/bin/time -l "$@" >"${log}.stdout" 2>"${log}" ;;
    *) "$@" ;;
  esac
}

giab_parse_time_log() {
  # Args: time_log_path → prints JSON fields wall_sec, max_rss_kb to stdout via python
  local log="$1"
  python3 - "$log" <<'PY'
import json, re, sys
text = open(sys.argv[1], encoding="utf-8", errors="replace").read()
wall = None
rss_kb = None
# GNU time -v
m = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*(\S+)", text)
if m:
    parts = m.group(1).split(":")
    if len(parts) == 3:
        wall = int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
    elif len(parts) == 2:
        wall = int(parts[0]) * 60 + float(parts[1])
m = re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", text)
if m:
    rss_kb = int(m.group(1))
# macOS /usr/bin/time -l formats (Darwin):
#   "        0.29 real         0.00 user         0.00 sys"
#   "             9715712  maximum resident set size"  (bytes, number BEFORE label)
# Older docs sometimes show "real 1.23" / "maximum resident set size  12345678".
if wall is None:
    # Darwin may use a locale decimal comma ("0,01 real").
    m = re.search(r"(\d+[.,]\d+)\s+real\b", text) or re.search(
        r"^\s*real\s+(\d+[.,]\d+)", text, re.M
    )
    if m:
        wall = float(m.group(1).replace(",", "."))
if rss_kb is None:
    m = re.search(r"(\d+)\s+maximum resident set size\b", text) or re.search(
        r"maximum resident set size\s+(\d+)", text
    )
    if m:
        # macOS reports bytes
        rss_kb = int(m.group(1)) // 1024
print(json.dumps({"wall_sec": wall, "max_rss_kb": rss_kb}))
PY
}

# hs37d5 autosomal lengths (numeric contigs; matches parity/realworld/assets/hs37d5.simple.fa).
# Used so ci-subset 50kb probes never exceed contig ends (GATK USER ERROR otherwise).
giab_hs37d5_chrom_len() {
  case "$1" in
    1) echo 249250621 ;;
    2) echo 243199373 ;;
    3) echo 198022430 ;;
    4) echo 191154276 ;;
    5) echo 180915260 ;;
    6) echo 171115067 ;;
    7) echo 159138663 ;;
    8) echo 146364022 ;;
    9) echo 141213431 ;;
    10) echo 135534747 ;;
    11) echo 135006516 ;;
    12) echo 133851895 ;;
    13) echo 115169878 ;;
    14) echo 107349540 ;;
    15) echo 102531392 ;;
    16) echo 90354753 ;;
    17) echo 81195210 ;;
    18) echo 78077248 ;;
    19) echo 59128983 ;;
    20) echo 63025520 ;;
    21) echo 48129895 ;;
    22) echo 51304566 ;;
    *)
      echo "[giab] unknown contig for hs37d5 length: $1" >&2
      return 2
      ;;
  esac
}

# Emit a 50kb probe for $1 (chrom), preferring start=chr*3e6+1e7, clamped into contig.
giab_ci_subset_probe() {
  local chr="$1"
  local probe_len=50000
  local len sample_start max_start end
  len="$(giab_hs37d5_chrom_len "${chr}")" || return 2
  sample_start=$((chr * 3000000 + 10000000))
  max_start=$((len - probe_len + 1))
  if (( max_start < 1 )); then
    echo "[giab] contig ${chr} shorter than ${probe_len} bp" >&2
    return 2
  fi
  if (( sample_start > max_start )); then
    # chr19/22 overflow the naive formula on hs37d5 — park at contig end.
    sample_start="${max_start}"
  fi
  if (( sample_start < 1 )); then
    sample_start=1
  fi
  end=$((sample_start + probe_len - 1))
  echo "${chr}:${sample_start}-${end}"
}

giab_build_intervals() {
  # Args: mode → writes interval strings one per line on stdout
  local mode="$1"
  case "${mode}" in
    smoke)
      # M4 / PR-smoke: small windows only
      echo "20:10000000-10050000"
      echo "21:41200001-41250000"
      echo "2:92300000-92350000"
      ;;
    ci-subset)
      # Practical “genome-wide” for CI: full chr20+chr21 + 50kb probes on other autosomes.
      # Probes are clamped to hs37d5 contig lengths (naive mid-ish formula overflows 19/22).
      echo "20:1-63025520"
      echo "21:1-48129895"
      local chr
      for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 22; do
        giab_ci_subset_probe "${chr}"
      done
      ;;
    chr20-21)
      echo "20:1-63025520"
      echo "21:1-48129895"
      ;;
    autosomes)
      local chr
      for chr in $(seq 1 22); do
        echo "${chr}"
      done
      ;;
    *)
      echo "[giab] unknown GIAB_MODE=${mode}" >&2
      return 2
      ;;
  esac
}

giab_mode_description() {
  local mode="$1"
  case "${mode}" in
    smoke)
      echo "SMOKE: three ~50kb windows (chr20/chr21/P12). Not genome-wide."
      ;;
    ci-subset)
      echo "CI-SUBSET (default “genome-wide” in this repo): FULL chr20 + FULL chr21 + one 50kb probe on each other autosome. Not all bases of chr1–19/22."
      ;;
    chr20-21)
      echo "CHR20-21: full chromosomes 20 and 21 only."
      ;;
    autosomes)
      echo "AUTOSOMES: chromosomes 1–22 in full. Requires large RAM/disk/time — not for M4 laptops."
      ;;
  esac
}
