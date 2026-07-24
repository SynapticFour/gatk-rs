#!/usr/bin/env bash
# Stage a *dense* NA12878 b37 BAM slice that overlaps GIAB high-confidence regions
# (Rust-native R3 — path out of L6 vacuous F1 / W-L6-1 on the sparse 20k BAM).
#
# Default source: NIST GIAB HG001 HiSeq 30× downsample (remote HTTPS + local BAI).
# Default window: chr20:10000000-10050000 (known GIAB-rich; small enough for CI/local).
#
# Usage:
#   ./scripts/parity/realworld/04_stage_na12878_giab_window_bam.sh
#
# Env:
#   J6_DENSE_INTERVAL   default 20:10000000-10050000
#   J6_DENSE_BAM_URL    remote BAM (Accept-Ranges required)
#   J6_DENSE_BAI_URL    remote/local BAI URL (downloaded once)
#   J6_DENSE_OUT_DIR    output directory
#
# Does not replace NA12878_20k (P12 spine parity). Wire into L6 via:
#   export P12_BAM=.../NA12878_giab_window.b37.bam
#   export P12_INTERVAL="$J6_DENSE_INTERVAL"
#   export J6_INTERVAL="$J6_DENSE_INTERVAL"
#   export J6_PARITY_INTERVAL="$J6_DENSE_INTERVAL"   # optional: skip spine lock
#   ./scripts/parity/run_hc_full_parity_j6_truth.sh
set -euo pipefail

export PATH="/usr/bin:/bin:/usr/sbin:/sbin:/opt/homebrew/bin:${PATH:-}"

script_dir="$(cd "$(dirname "$0")" && pwd)"
repo_root="$(cd "${script_dir}/../../.." && pwd)"
cd "${repo_root}"

interval="${J6_DENSE_INTERVAL:-20:10000000-10050000}"
out_dir="${J6_DENSE_OUT_DIR:-${repo_root}/parity/realworld/na12878_giab_window_b37}"
bam_url="${J6_DENSE_BAM_URL:-https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam}"
bai_url="${J6_DENSE_BAI_URL:-${bam_url}.bai}"

mkdir -p "${out_dir}"
local_bai="${out_dir}/RMNISTHS_30xdownsample.bam.bai"
out_bam="${out_dir}/NA12878_giab_window.b37.bam"
out_bai="${out_dir}/NA12878_giab_window.b37.bai"
meta="${out_dir}/stage_meta.json"

if ! command -v samtools >/dev/null 2>&1; then
  echo "[j6-dense] samtools required on PATH" >&2
  exit 1
fi
if ! command -v curl >/dev/null 2>&1; then
  echo "[j6-dense] curl required on PATH" >&2
  exit 1
fi

echo "=== 04_stage_na12878_giab_window_bam ==="
echo "interval=${interval}"
echo "out_bam=${out_bam}"

if [[ ! -f "${local_bai}" ]]; then
  echo "[j6-dense] downloading BAI ($(basename "${local_bai}"))…"
  curl -fsSL "${bai_url}" -o "${local_bai}.partial"
  mv -f "${local_bai}.partial" "${local_bai}"
fi

need_slice=1
if [[ -f "${out_bam}" && -f "${out_bai}" && -f "${meta}" ]]; then
  staged_iv="$(python3 -c "import json; print(json.load(open('${meta}')).get('interval',''))" 2>/dev/null || true)"
  if [[ "${staged_iv}" == "${interval}" ]]; then
    need_slice=0
    echo "[j6-dense] reuse existing slice for ${interval}"
  fi
fi

if [[ "${need_slice}" == "1" ]]; then
  echo "[j6-dense] slicing remote BAM (HTTPS range via samtools)…"
  # -X: BAM URL + local BAI
  samtools view -b -X "${bam_url}" "${local_bai}" "${interval}" > "${out_bam}.partial"
  mv -f "${out_bam}.partial" "${out_bam}"
  samtools index "${out_bam}" "${out_bai}"
  python3 - "${meta}" "${interval}" "${bam_url}" "${out_bam}" <<'PY'
import json, pathlib, sys, subprocess
meta, interval, url, bam = sys.argv[1:5]
n = int(subprocess.check_output(["samtools", "view", "-c", bam], text=True).strip())
pathlib.Path(meta).write_text(
    json.dumps(
        {
            "label": "na12878-giab-window-b37",
            "interval": interval,
            "source_bam_url": url,
            "bam": bam,
            "read_count": n,
            "notes": "Dense slice for non-vacuous GIAB truth (R3). Not the P12 spine corpus.",
        },
        indent=2,
    )
    + "\n",
    encoding="utf-8",
)
print(f"[j6-dense] reads_in_slice={n}")
PY
fi

reads="$(samtools view -c "${out_bam}")"
echo "[j6-dense] staged BAM read_count=${reads}"
if [[ "${reads}" -lt 100 ]]; then
  echo "[j6-dense] WARNING: slice looks sparse (reads<100); try a different J6_DENSE_INTERVAL" >&2
fi

echo "export P12_BAM=\"${out_bam}\""
echo "export P12_INTERVAL=\"${interval}\""
echo "export J6_INTERVAL=\"${interval}\""
echo "04_stage_na12878_giab_window_bam: done"
