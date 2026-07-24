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
      # Practical “genome-wide” for CI: full chr20+chr21 + 50kb samples on other autosomes.
      echo "20:1-63025520"
      echo "21:1-48129895"
      local chr sample_start
      for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 22; do
        # Deterministic mid-chromosome-ish 50kb probes (hs37d5 lengths not required for -L).
        sample_start=$((chr * 3000000 + 10000000))
        echo "${chr}:${sample_start}-$((sample_start + 49999))"
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
