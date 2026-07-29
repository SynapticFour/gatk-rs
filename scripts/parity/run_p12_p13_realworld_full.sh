#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
report_dir="${repo_root}/parity/reports"
data_dir="${P12P13_DATA_DIR:-${repo_root}/parity/realworld/assets}"
mkdir -p "${report_dir}" "${data_dir}"
cd "${repo_root}"

log_file="${report_dir}/p12_p13_realworld_full.log"
timestamp="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"

log() {
  local msg="$1"
  printf '[%s] %s\n' "$(date -u +"%Y-%m-%dT%H:%M:%SZ")" "${msg}" | tee -a "${log_file}"
}

run_cmd() {
  log "RUN: $*"
  "$@" 2>&1 | tee -a "${log_file}"
}

log "=== P12/P13 real-world full run start (${timestamp}) ==="

# Prefer NCBI mirrors from GitHub-hosted runners; EBI often stalls for hours.
# Best path: GitHub Release CDN (giab-ref-v1) — no FTP at all.
ref_gz_url="${P12P13_REF_GZ_URL:-}"
ref_gz="${data_dir}/hs37d5.fa.gz"
ref_fa_raw="${data_dir}/hs37d5.fa"
ref_fa="${data_dir}/hs37d5.simple.fa"
ref_fai="${ref_fa}.fai"
ref_dict="${data_dir}/hs37d5.simple.dict"

truth_url="${P12P13_TRUTH_URL:-https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh37/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz}"
truth_vcf="${data_dir}/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
truth_bed_url="${P12P13_TRUTH_BED_URL:-https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh37/HG001_GRCh37_1_22_v4.2.1_benchmark.bed}"
truth_bed="${data_dir}/HG001_GRCh37_1_22_v4.2.1_benchmark.bed"

# --- Preferred: pinned release assets (simple.fa.gz + fai + dict) -------------
if [[ ! -f "${ref_fa}" || ! -f "${ref_fai}" || ! -f "${ref_dict}" ]]; then
  if [[ "${GIAB_REF_USE_RELEASE:-1}" == "1" ]]; then
    log "Trying GitHub Release giab-ref assets…"
    if "${repo_root}/scripts/parity/giab/fetch_hs37d5_release.sh"; then
      log "REF from GitHub Release OK"
    else
      log "GitHub Release fetch failed; will try FTP mirrors"
    fi
  fi
fi

download_ref_gz() {
  local mirrors=()
  if [[ -n "${ref_gz_url}" ]]; then
    mirrors+=("${ref_gz_url}")
  fi
  mirrors+=(
    "https://ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
    "https://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
  )
  local url attempt
  rm -f "${ref_gz}.partial"
  for url in "${mirrors[@]}"; do
    log "REF_GZ trying ${url}"
    for attempt in 1 2 3; do
      # Resume partial downloads; NCBI is typically reachable from GH runners.
      if curl -fL --retry 2 --retry-delay 5 --connect-timeout 20 --max-time 1800 \
        -C - "${url}" -o "${ref_gz}.partial"; then
        # Require plausible size (~850–950 MiB) + valid gzip — small HTML error
        # bodies used to pass curl -f and then fail later as "corrupt gzip".
        local sz
        sz="$(wc -c < "${ref_gz}.partial" | tr -d ' ')"
        if [[ "${sz}" -lt 500000000 ]]; then
          log "REF_GZ too small (${sz} bytes) from ${url}; retrying"
          rm -f "${ref_gz}.partial"
        elif gzip -t "${ref_gz}.partial" 2>/dev/null; then
          mv "${ref_gz}.partial" "${ref_gz}"
          log "REF_GZ downloaded from ${url} (${sz} bytes)"
          return 0
        else
          log "REF_GZ corrupt gzip from ${url}; retrying"
          rm -f "${ref_gz}.partial"
        fi
      else
        log "REF_GZ curl failed (${url} attempt ${attempt})"
      fi
      sleep $((attempt * 5))
    done
  done
  return 1
}

# --- Fallback FTP path (only if release did not provide simple.fa) -------------
if [[ ! -f "${ref_fa}" ]]; then
  if [[ ! -f "${ref_gz}" ]]; then
    if ! download_ref_gz; then
      log "REF_GZ download failed on all mirrors"
      exit 1
    fi
  else
    log "REF_GZ already present: ${ref_gz}"
  fi

  if [[ ! -f "${ref_fa_raw}" ]]; then
    run_cmd bash -lc "gzip -dc '${ref_gz}' > '${ref_fa_raw}'"
  else
    log "REF_FASTA_RAW already present: ${ref_fa_raw}"
  fi

  if [[ ! -f "${ref_fa}" ]]; then
    run_cmd python3 - "${ref_fa_raw}" "${ref_fa}" <<'PY'
import pathlib
import sys
src = pathlib.Path(sys.argv[1])
dst = pathlib.Path(sys.argv[2])
with src.open("r", encoding="utf-8", errors="replace") as fin, dst.open("w", encoding="utf-8") as fout:
    for line in fin:
        if line.startswith(">"):
            token = line[1:].strip().split()[0]
            fout.write(f">{token}\n")
        else:
            fout.write(line)
print(dst)
PY
  else
    log "REF_FASTA_SIMPLE already present: ${ref_fa}"
  fi
else
  log "REF_FASTA_SIMPLE already present: ${ref_fa}"
fi

if [[ ! -f "${ref_fai}" ]]; then
  # Prefer host samtools (CI installs it) — avoids a multi‑GB GATK Docker pull
  # just to faidx the reference during GIAB prepare.
  if command -v samtools >/dev/null 2>&1; then
    run_cmd samtools faidx "${ref_fa}"
  else
    run_cmd docker run --rm --platform "${GATK_DOCKER_PLATFORM:-linux/amd64}" \
      -v "${repo_root}:${repo_root}" \
      -w "${repo_root}" \
      "${GATK_DOCKER_IMAGE:-us.gcr.io/broad-gatk/gatk:4.4.0.0}" \
      samtools faidx "${ref_fa}"
  fi
else
  log "REF_FAI already present: ${ref_fai}"
fi

if [[ ! -f "${ref_dict}" ]]; then
  if command -v samtools >/dev/null 2>&1; then
    run_cmd samtools dict -o "${ref_dict}" "${ref_fa}"
  else
    run_cmd docker run --rm --platform "${GATK_DOCKER_PLATFORM:-linux/amd64}" \
      -v "${repo_root}:${repo_root}" \
      -w "${repo_root}" \
      "${GATK_DOCKER_IMAGE:-us.gcr.io/broad-gatk/gatk:4.4.0.0}" \
      gatk CreateSequenceDictionary -R "${ref_fa}" -O "${ref_dict}" --QUIET true
  fi
else
  log "REF_DICT already present: ${ref_dict}"
fi

if [[ ! -f "${truth_vcf}" ]]; then
  run_cmd curl -fsSL "${truth_url}" -o "${truth_vcf}"
else
  log "TRUTH already present: ${truth_vcf}"
fi

if [[ ! -f "${truth_bed}" ]]; then
  run_cmd curl -fsSL "${truth_bed_url}" -o "${truth_bed}"
else
  log "TRUTH_BED already present: ${truth_bed}"
fi

# Playbook: fetch hs37d5 + GIAB only, then exit (no BAM/HC). See docs/realworld-parity-playbook.md
if [[ "${REALWORLD_STOP_AFTER_ASSETS:-0}" == "1" ]]; then
  log "REALWORLD_STOP_AFTER_ASSETS=1 — stopping after reference + truth staging."
  exit 0
fi

if [[ -z "${P12_INTERVAL:-}" ]]; then
  # Choose a window with observed reads so runtime stays bounded but meaningful.
  chrom="${P12_CHROM:-20}"
  win="${P12_WINDOW_SIZE:-50000}"
  first_pos="$(
    docker run --rm --platform "${GATK_DOCKER_PLATFORM:-linux/amd64}" \
      -v "${repo_root}:${repo_root}" \
      -w "${repo_root}" \
      "${GATK_DOCKER_IMAGE:-us.gcr.io/broad-gatk/gatk:4.4.0.0}" \
      samtools view "${repo_root}/parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam" "${chrom}" \
      | awk 'NR==1 {print $4; exit}' || true
  )"
  if [[ -z "${first_pos}" ]]; then
    first_pos=1
  fi
  start="${first_pos}"
  end="$((start + win - 1))"
  export P12_INTERVAL="${chrom}:${start}-${end}"
  log "AUTO_INTERVAL selected from BAM coverage: ${P12_INTERVAL}"
else
  export P12_INTERVAL="${P12_INTERVAL}"
fi

export P12_REFERENCE="${ref_fa}"
export P13_TRUTH_VCF="${truth_vcf}"
export P13_CHROM="${P13_CHROM:-20}"
export P13_REGIONS_BED="${P13_REGIONS_BED:-${truth_bed}}"

run_cmd ./scripts/parity/run_p12_realworld_na12878_20k.sh
run_cmd ./scripts/parity/run_p13_realworld_truth_eval.sh

summary_json="${report_dir}/p12_p13_realworld_full_summary.json"
summary_md="${report_dir}/p12_p13_realworld_full_summary.md"
run_cmd python3 - "${report_dir}" "${summary_json}" "${summary_md}" "${log_file}" <<'PY'
import json
import pathlib
import sys

reports = pathlib.Path(sys.argv[1])
summary_json = pathlib.Path(sys.argv[2])
summary_md = pathlib.Path(sys.argv[3])
log_file = pathlib.Path(sys.argv[4])

p12 = json.loads((reports / "p12_realworld_na12878_20k.json").read_text(encoding="utf-8"))
p13 = json.loads((reports / "p13_realworld_truth_eval.json").read_text(encoding="utf-8"))

payload = {
    "label": "phase12-13-realworld-full-run",
    "status": "pass" if p12.get("status") == "pass" and p13.get("status") == "pass" else "needs_attention",
    "p12_status": p12.get("status"),
    "p13_status": p13.get("status"),
    "p12": p12,
    "p13": p13,
    "log_file": str(log_file),
}
summary_json.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
summary_md.write_text(
    "\n".join(
        [
            "# P12/P13 Real-world Full Run Summary",
            "",
            f"- p12 status: **{p12.get('status')}**",
            f"- p12 parity_status: **{p12.get('parity_status', 'unknown')}**",
            f"- p13 status: **{p13.get('status')}**",
            f"- p13 eval_interval: `{p13.get('eval_interval')}`",
            f"- p12 java/rust variants: `{p12.get('java_variant_count')} / {p12.get('rust_variant_count')}`",
            f"- p13 java F1: `{(p13.get('java') or {}).get('f1', 0):.6f}`",
            f"- p13 rust F1: `{(p13.get('rust') or {}).get('f1', 0):.6f}`",
            f"- run log: `{log_file}`",
        ]
    )
    + "\n",
    encoding="utf-8",
)
print(f"[p12p13-full] wrote {summary_json}")
print(f"[p12p13-full] wrote {summary_md}")
PY

log "=== P12/P13 real-world full run complete ==="
