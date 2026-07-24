#!/usr/bin/env bash
# Download official GIAB truth VCFs/BEDs (HG001/HG002/HG005) + key GRCh37 stratification BEDs.
# Does NOT download multi-GB stratification tarballs or WGS BAMs (see run_genomewide_equivalence.sh).
#
# Usage:
#   ./scripts/parity/giab/fetch_giab_truthsets.sh
#
# Env:
#   GIAB_TRUTH_ROOT   default: parity/giab/truth
#   GIAB_STRAT_ROOT   default: parity/giab/stratifications/GRCh37
#   GIAB_FORCE=1      re-download even if present
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../../.." && pwd)"
# shellcheck source=../m4_disk_guard.sh
source "${repo_root}/scripts/parity/m4_disk_guard.sh"
m4_require_free_gb "${GIAB_FETCH_MIN_FREE_GB:-6}" || exit 1

truth_root="${GIAB_TRUTH_ROOT:-${repo_root}/parity/giab/truth}"
strat_root="${GIAB_STRAT_ROOT:-${repo_root}/parity/giab/stratifications/GRCh37}"
ftp_base="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab"
strat_base="${ftp_base}/release/genome-stratifications/v3.3/GRCh37@all"
force="${GIAB_FORCE:-0}"

mkdir -p "${truth_root}" "${strat_root}"

download() {
  local url="$1" dest="$2"
  if [[ -f "${dest}" && "${force}" != "1" ]]; then
    echo "[giab-fetch] reuse $(basename "${dest}")"
    return 0
  fi
  echo "[giab-fetch] GET ${url}"
  curl -fL --retry 3 --retry-delay 2 -o "${dest}.partial" "${url}"
  mv -f "${dest}.partial" "${dest}"
}

echo "=== fetch_giab_truthsets ==="
echo "truth_root=${truth_root}"
echo "strat_root=${strat_root}"

# --- Truth callsets (NISTv4.2.1 / GRCh37) ------------------------------------
# HG001
download \
  "${ftp_base}/release/NA12878_HG001/NISTv4.2.1/GRCh37/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz" \
  "${truth_root}/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
download \
  "${ftp_base}/release/NA12878_HG001/NISTv4.2.1/GRCh37/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz.tbi" \
  "${truth_root}/HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz.tbi" || true
download \
  "${ftp_base}/release/NA12878_HG001/NISTv4.2.1/GRCh37/HG001_GRCh37_1_22_v4.2.1_benchmark.bed" \
  "${truth_root}/HG001_GRCh37_1_22_v4.2.1_benchmark.bed"

# HG002 (confident regions file name differs: *_noinconsistent.bed)
download \
  "${ftp_base}/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz" \
  "${truth_root}/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
download \
  "${ftp_base}/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz.tbi" \
  "${truth_root}/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz.tbi" || true
download \
  "${ftp_base}/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh37/HG002_GRCh37_1_22_v4.2.1_benchmark_noinconsistent.bed" \
  "${truth_root}/HG002_GRCh37_1_22_v4.2.1_benchmark.bed"

# HG005
download \
  "${ftp_base}/release/ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh37/HG005_GRCh37_1_22_v4.2.1_benchmark.vcf.gz" \
  "${truth_root}/HG005_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
download \
  "${ftp_base}/release/ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh37/HG005_GRCh37_1_22_v4.2.1_benchmark.vcf.gz.tbi" \
  "${truth_root}/HG005_GRCh37_1_22_v4.2.1_benchmark.vcf.gz.tbi" || true
download \
  "${ftp_base}/release/ChineseTrio/HG005_NA24631_son/NISTv4.2.1/GRCh37/HG005_GRCh37_1_22_v4.2.1_benchmark.bed" \
  "${truth_root}/HG005_GRCh37_1_22_v4.2.1_benchmark.bed"

# --- Stratifications (individual BEDs — not the ~1.4GiB @all tarball) --------
# Source: https://github.com/genome-in-a-bottle/genome-stratifications (FTP mirror v3.3)
strat_files=(
  "LowComplexity/GRCh37_AllTandemRepeatsandHomopolymers_slop5.bed.gz"
  "SegmentalDuplications/GRCh37_segdups.bed.gz"
  "Union/GRCh37_alldifficultregions.bed.gz"
  "Union/GRCh37_notinalldifficultregions.bed.gz"
  "OtherDifficult/GRCh37_MHC.bed.gz"
)

for rel in "${strat_files[@]}"; do
  dest="${strat_root}/$(basename "${rel}")"
  download "${strat_base}/${rel}" "${dest}"
done

# Manifest for consumers
python3 - "${truth_root}" "${strat_root}" <<'PY'
import json, pathlib, sys
truth, strat = map(pathlib.Path, sys.argv[1:3])
manifest = {
    "assembly": "GRCh37",
    "benchmark": "NISTv4.2.1",
    "samples": {
        "HG001": {
            "truth_vcf": str(truth / "HG001_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"),
            "confident_bed": str(truth / "HG001_GRCh37_1_22_v4.2.1_benchmark.bed"),
        },
        "HG002": {
            "truth_vcf": str(truth / "HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"),
            "confident_bed": str(truth / "HG002_GRCh37_1_22_v4.2.1_benchmark.bed"),
            "notes": "BED is official *_noinconsistent.bed renamed for uniform layout",
        },
        "HG005": {
            "truth_vcf": str(truth / "HG005_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"),
            "confident_bed": str(truth / "HG005_GRCh37_1_22_v4.2.1_benchmark.bed"),
        },
    },
    "stratifications": sorted(str(p) for p in strat.glob("*.bed.gz")),
    "stratification_source": "ftp genome-stratifications v3.3 GRCh37@all (selected beds)",
}
(truth / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
print(f"[giab-fetch] wrote {truth / 'manifest.json'}")
PY

echo "[giab-fetch] done"
