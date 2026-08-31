#!/usr/bin/env bash
# 6R.43 — freeze inputs and run pinned GATK 4.4.0.0 + gatk-rs HC on the holdout panel.
# Does not modify production algorithms.
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
# shellcheck source=docs/GATK_PINNED.env
source "${repo_root}/docs/GATK_PINNED.env"

panel="${repo_root}/scripts/parity/6r43_holdout_panel.json"
out_root="${repo_root}/parity/reports/6r43"
gatk_image="${GATK_DOCKER_IMAGE:-broadinstitute/gatk:4.4.0.0}"
gatk_platform="${GATK_DOCKER_PLATFORM:-linux/amd64}"
rust_rev="$(git rev-parse HEAD)"
only="${HOLDOUT_ONLY:-}"

mkdir -p "${out_root}"
echo "${rust_rev}" > "${out_root}/RUST_REV"
echo "${GATK_PINNED_SHA}" > "${out_root}/JAVA_SHA"

ref="${repo_root}/parity/realworld/assets/hs37d5.simple.fa"
if [[ ! -f "${ref}.dict" && -f "${repo_root}/parity/realworld/assets/hs37d5.simple.dict" ]]; then
  ln -sfn "${repo_root}/parity/realworld/assets/hs37d5.simple.dict" "${ref}.dict"
fi

run_java() {
  local id="$1" interval="$2" bam="$3" dest="$4"
  mkdir -p "$(dirname "${dest}")"
  if [[ -f "${dest}" && "${HOLDOUT_FORCE:-0}" != "1" ]]; then
    echo "[6R.43] reuse java ${id}"
    return 0
  fi
  echo "[6R.43] JAVA ${id} ${interval}"
  docker run --rm --platform "${gatk_platform}" \
    -v "${repo_root}:${repo_root}" \
    -w "${repo_root}" \
    "${gatk_image}" \
    gatk --java-options "-Xmx4g" HaplotypeCaller \
    -R "${ref}" \
    -I "${repo_root}/${bam}" \
    -O "${dest}" \
    -L "${interval}" \
    --native-pair-hmm-threads 1 \
    --verbosity ERROR
}

ids="$(python3 - "${panel}" "${only}" <<'PY'
import json, sys
panel = json.load(open(sys.argv[1]))
only = sys.argv[2]
for r in panel["regions"]:
    if only and r["id"] != only:
        continue
    print("\t".join([r["id"], r["interval"], r["bam"]]))
PY
)"

while IFS=$'\t' read -r id interval bam; do
  [[ -z "${id}" ]] && continue
  bam_abs="${repo_root}/${bam}"
  if [[ ! -f "${bam_abs}" ]]; then
    echo "[6R.43] SKIP ${id}: missing BAM ${bam}" >&2
    continue
  fi
  dest="${out_root}/${id}/java.vcf"
  run_java "${id}" "${interval}" "${bam}" "${dest}"
done <<<"${ids}"

echo "[6R.43] java pass complete → ${out_root}"
