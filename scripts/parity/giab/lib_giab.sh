#!/usr/bin/env bash
# Shared helpers for GIAB genome-wide equivalence runners.
# shellcheck shell=bash

# Resolved at source time (inside functions BASH_SOURCE[0] is the caller).
_GIAB_LIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

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

# P12 spine interval used by GIAB smoke (and P12 L* gates).
GIAB_SMOKE_P12_INTERVAL="2:92300000-92350000"

# Public NA12878_20k_b37 corpus (same evidence class as P12 L* gates).
GIAB_P12_20K_S3_BASE_URL="${GIAB_P12_20K_S3_BASE_URL:-https://gatk-test-data.s3.amazonaws.com/wgs_bam/NA12878_20k_b37}"

giab_ensure_na12878_20k_bam() {
  # Download NA12878_20k_b37 into parity/realworld if missing. Prints BAM path on stdout.
  local repo_root="$1"
  local data_dir="${GIAB_P12_20K_DIR:-${repo_root}/parity/realworld/na12878_20k_b37}"
  local bam="${data_dir}/NA12878_20k.b37.bam"
  local bai="${data_dir}/NA12878_20k.b37.bai"
  mkdir -p "${data_dir}"
  if [[ ! -f "${bam}" ]]; then
    echo "[giab] download NA12878_20k.b37.bam (smoke P12 evidence class)…" >&2
    curl -fL --retry 3 -o "${bam}.partial" "${GIAB_P12_20K_S3_BASE_URL}/NA12878_20k.b37.bam"
    mv -f "${bam}.partial" "${bam}"
  fi
  if [[ ! -f "${bai}" ]]; then
    curl -fL --retry 3 -o "${bai}.partial" "${GIAB_P12_20K_S3_BASE_URL}/NA12878_20k.b37.bai"
    mv -f "${bai}.partial" "${bai}"
  fi
  printf '%s\n' "${bam}"
}

giab_stage_smoke_bam_hybrid() {
  # Stage smoke BAM: chr20/chr21 from full HG001 30× URL; P12 spine from NA12878_20k.
  # Full-30× P12 is centromere-scale (~537k reads after positional DS) and is a
  # benchmark-host gate — not safe on 16 GiB hosted runners.
  # Args: sample bam_url out_bam repo_root intervals...
  local sample="$1" bam_url="$2" out_bam="$3" repo_root="$4"
  shift 4
  local -a ivs=("$@")
  local bai_local="${out_bam}.remote.bai"
  local out_bai="${out_bam}.bai"
  local tmp_dir normal_bam p12_bam p12_src iv
  local -a normal_ivs=()
  local have_p12=0
  mkdir -p "$(dirname "${out_bam}")"
  if [[ -f "${out_bam}" && -f "${out_bai}" ]]; then
    echo "[giab] reuse smoke hybrid BAM ${out_bam}"
    return 0
  fi
  for iv in "${ivs[@]}"; do
    if [[ "${iv}" == "${GIAB_SMOKE_P12_INTERVAL}" ]]; then
      have_p12=1
    else
      normal_ivs+=("${iv}")
    fi
  done
  tmp_dir="$(mktemp -d "${TMPDIR:-/tmp}/giab-smoke-hybrid.XXXXXX")"
  normal_bam="${tmp_dir}/normal.bam"
  p12_bam="${tmp_dir}/p12.bam"
  if [[ "${#normal_ivs[@]}" -gt 0 ]]; then
    if [[ ! -f "${bai_local}" ]]; then
      echo "[giab] download BAI for ${sample}…"
      curl -fL --retry 3 -o "${bai_local}.partial" "${bam_url}.bai"
      mv -f "${bai_local}.partial" "${bai_local}"
    fi
    echo "[giab] slicing remote BAM for ${sample} smoke non-P12 (${#normal_ivs[@]} intervals)…"
    samtools view -b -X "${bam_url}" "${bai_local}" "${normal_ivs[@]}" > "${normal_bam}"
  fi
  if [[ "${have_p12}" -eq 1 ]]; then
    p12_src="$(giab_ensure_na12878_20k_bam "${repo_root}")"
    echo "[giab] slicing NA12878_20k for smoke P12 (${GIAB_SMOKE_P12_INTERVAL})…"
    samtools view -b "${p12_src}" "${GIAB_SMOKE_P12_INTERVAL}" > "${p12_bam}"
  fi
  if [[ "${#normal_ivs[@]}" -gt 0 && "${have_p12}" -eq 1 ]]; then
    samtools merge -f "${out_bam}.partial" "${normal_bam}" "${p12_bam}"
  elif [[ "${#normal_ivs[@]}" -gt 0 ]]; then
    mv -f "${normal_bam}" "${out_bam}.partial"
  elif [[ "${have_p12}" -eq 1 ]]; then
    mv -f "${p12_bam}" "${out_bam}.partial"
  else
    rm -rf "${tmp_dir}"
    echo "[giab] stage_smoke_bam_hybrid: empty interval list" >&2
    return 2
  fi
  mv -f "${out_bam}.partial" "${out_bam}"
  samtools index "${out_bam}" "${out_bai}"
  rm -rf "${tmp_dir}"
}

giab_build_intervals() {
  # Args: mode → writes interval strings one per line on stdout
  local mode="$1"
  case "${mode}" in
    smoke)
      # M4 / PR-smoke: small windows only (P12 spine staged from NA12878_20k; see hybrid BAM)
      echo "20:10000000-10050000"
      echo "21:41200001-41250000"
      echo "${GIAB_SMOKE_P12_INTERVAL}"
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
      echo "SMOKE: three ~50kb windows (chr20/chr21/P12). P12 reads from NA12878_20k evidence class; chr20/21 from HG001 30×. Full-30× P12 is benchmark-host only. Not genome-wide."
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

# Window size (bp) for splitting full-contig intervals into matrix jobs.
# Contig-sized shards (full chr20/21) can still burn the GitHub-hosted 6h hard
# cap at GIAB_THREADS=2; 10 Mb windows keep each shard×engine job smaller.
giab_hc_window_bp() {
  echo "${GIAB_HC_WINDOW_BP:-10000000}"
}

# Parse "chr", "chr:start-end" → chrom start end (1-based inclusive).
giab_parse_interval() {
  local iv="$1"
  local chrom rest start end len
  case "${iv}" in
    *:*)
      chrom="${iv%%:*}"
      rest="${iv#*:}"
      start="${rest%-*}"
      end="${rest#*-}"
      ;;
    *)
      chrom="${iv}"
      len="$(giab_hs37d5_chrom_len "${chrom}")" || return 2
      start=1
      end="${len}"
      ;;
  esac
  printf '%s %s %s\n' "${chrom}" "${start}" "${end}"
}

# Write one .intervals file per window for a full-chrom (or long) interval.
# Names: ${prefix}_w00, ${prefix}_w01, … (zero-padded for concat order).
# Short intervals (≤ window_bp) become a single ${prefix}_w00 shard.
giab_write_interval_windows() {
  local shard_dir="$1"
  local prefix="$2"
  local iv="$3"
  local window_bp="$4"
  local chrom start end wstart wend idx
  read -r chrom start end < <(giab_parse_interval "${iv}") || return 2
  if (( start < 1 || end < start )); then
    echo "[giab] bad interval ${iv}" >&2
    return 2
  fi
  idx=0
  wstart="${start}"
  while (( wstart <= end )); do
    wend=$((wstart + window_bp - 1))
    if (( wend > end )); then
      wend="${end}"
    fi
    printf '%s:%s-%s\n' "${chrom}" "${wstart}" "${wend}" \
      > "${shard_dir}/$(printf '%s_w%02d' "${prefix}" "${idx}").intervals"
    wstart=$((wend + 1))
    idx=$((idx + 1))
  done
}

# Write shard interval lists under $1 from $2 (intervals.txt). One shard per file:
#   smoke        → 00_all
#   ci-subset    → 00_chr20_wNN, 01_chr21_wNN, 02_probes  (matrix jobs under 6h)
#   chr20-21     → 00_chr20_wNN, 01_chr21_wNN
#   autosomes    → one shard per contig (self-hosted; not window-split)
giab_write_hc_shards() {
  local shard_dir="$1"
  local intervals_file="$2"
  local mode="$3"
  local window_bp
  window_bp="$(giab_hc_window_bp)"
  mkdir -p "${shard_dir}"
  find "${shard_dir}" -maxdepth 1 -type f -name '*.intervals' -delete 2>/dev/null || true

  case "${mode}" in
    smoke)
      cp "${intervals_file}" "${shard_dir}/00_all.intervals"
      ;;
    ci-subset)
      : > "${shard_dir}/02_probes.intervals"
      while IFS= read -r iv || [[ -n "${iv}" ]]; do
        [[ -z "${iv}" ]] && continue
        case "${iv}" in
          20:* | 20)
            giab_write_interval_windows "${shard_dir}" "00_chr20" "${iv}" "${window_bp}"
            ;;
          21:* | 21)
            giab_write_interval_windows "${shard_dir}" "01_chr21" "${iv}" "${window_bp}"
            ;;
          *)
            echo "${iv}" >> "${shard_dir}/02_probes.intervals"
            ;;
        esac
      done < "${intervals_file}"
      ;;
    chr20-21)
      while IFS= read -r iv || [[ -n "${iv}" ]]; do
        [[ -z "${iv}" ]] && continue
        case "${iv}" in
          20:* | 20)
            giab_write_interval_windows "${shard_dir}" "00_chr20" "${iv}" "${window_bp}"
            ;;
          21:* | 21)
            giab_write_interval_windows "${shard_dir}" "01_chr21" "${iv}" "${window_bp}"
            ;;
          *)
            echo "[giab] unexpected interval for chr20-21 mode: ${iv}" >&2
            return 2
            ;;
        esac
      done < "${intervals_file}"
      ;;
    autosomes)
      local chr
      for chr in $(seq 1 22); do
        printf '%s\n' "${chr}" > "${shard_dir}/$(printf '%02d' "${chr}")_chr${chr}.intervals"
      done
      ;;
    *)
      echo "[giab] unknown mode for sharding: ${mode}" >&2
      return 2
      ;;
  esac

  # Drop empty shard files (e.g. missing probe file).
  find "${shard_dir}" -maxdepth 1 -type f -name '*.intervals' -size 0 -delete 2>/dev/null || true
}

# Concatenate VCF shards → $1. Prefers bcftools; falls back to concat_vcfs.py.
giab_concat_vcfs() {
  local out_vcf="$1"
  shift
  local -a inputs=("$@")
  if [[ "${#inputs[@]}" -eq 0 ]]; then
    echo "[giab] giab_concat_vcfs: no inputs" >&2
    return 2
  fi
  if [[ "${#inputs[@]}" -eq 1 ]]; then
    cp -f "${inputs[0]}" "${out_vcf}"
    return 0
  fi
  if command -v bcftools >/dev/null 2>&1; then
    bcftools concat -a -O v -o "${out_vcf}" "${inputs[@]}"
    return 0
  fi
  python3 "${_GIAB_LIB_DIR}/concat_vcfs.py" -o "${out_vcf}" "${inputs[@]}"
}
