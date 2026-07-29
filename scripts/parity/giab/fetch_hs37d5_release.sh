#!/usr/bin/env bash
# Fetch pinned hs37d5.simple FASTA (+ fai/dict) from the giab-ref GitHub Release.
# Used by GIAB prepare so GitHub-hosted runners never curl EBI/NCBI FTP.
#
# Override:
#   GIAB_REF_RELEASE_TAG=giab-ref-v1
#   GIAB_REF_RELEASE_BASE=https://github.com/OWNER/REPO/releases/download/TAG
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
data_dir="${P12P13_DATA_DIR:-${repo_root}/parity/realworld/assets}"
mkdir -p "${data_dir}"

tag="${GIAB_REF_RELEASE_TAG:-giab-ref-v1}"
base="${GIAB_REF_RELEASE_BASE:-https://github.com/SynapticFour/gatk-rs/releases/download/${tag}}"

ref_fa="${data_dir}/hs37d5.simple.fa"
ref_gz="${data_dir}/hs37d5.simple.fa.gz"
ref_fai="${ref_fa}.fai"
ref_dict="${data_dir}/hs37d5.simple.dict"
sums="${data_dir}/SHA256SUMS.giab-ref"

download() {
  local url="$1" out="$2"
  echo "[giab-ref] GET ${url}"
  curl -fL --retry 5 --retry-delay 3 --connect-timeout 20 --max-time 900 \
    -o "${out}.partial" "${url}"
  mv "${out}.partial" "${out}"
}

if [[ -f "${ref_fa}" && -f "${ref_fai}" && -f "${ref_dict}" ]]; then
  echo "[giab-ref] already present: ${ref_fa}"
  exit 0
fi

download "${base}/SHA256SUMS" "${sums}"
download "${base}/hs37d5.simple.fa.gz" "${ref_gz}"
download "${base}/hs37d5.simple.fa.fai" "${ref_fai}"
download "${base}/hs37d5.simple.dict" "${ref_dict}"

echo "[giab-ref] verifying checksums…"
(
  cd "${data_dir}"
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum -c SHA256SUMS.giab-ref
  else
    shasum -a 256 -c SHA256SUMS.giab-ref
  fi
)

echo "[giab-ref] decompressing ${ref_gz} → ${ref_fa}…"
gzip -dc "${ref_gz}" > "${ref_fa}.partial"
mv "${ref_fa}.partial" "${ref_fa}"

if command -v samtools >/dev/null 2>&1; then
  n_fai="$(wc -l < "${ref_fai}" | tr -d ' ')"
  echo "[giab-ref] fai contigs=${n_fai}"
fi

ls -lh "${ref_fa}" "${ref_fai}" "${ref_dict}"
echo "[giab-ref] done"
