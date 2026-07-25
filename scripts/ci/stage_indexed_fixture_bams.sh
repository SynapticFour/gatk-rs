#!/usr/bin/env bash
# Build samtools-indexed BAMs under parity/build/sam-indexed-bam/ from fixture SAMs.
# parity/build/ is gitignored; integration tests (e.g. indel_threading_graph) expect these caches.
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"

if ! command -v samtools >/dev/null 2>&1; then
  echo "[stage-bams] ERROR: samtools not on PATH" >&2
  exit 1
fi

mkdir -p parity/build/sam-indexed-bam
for base in \
  p5_live_case_indel \
  p5_live_case_snp \
  p5_live_case_ref \
  p5_live_case_repeat \
  p5_live_case_lowcomplex \
  p11_java_positive \
  p11_field_diff_case \
  read_filter_slice \
  read_order
do
  sam="parity/fixtures/${base}.sam"
  out="parity/build/sam-indexed-bam/${base}.bam"
  if [[ ! -f "${sam}" ]]; then
    echo "[stage-bams] missing ${sam}; skip"
    continue
  fi
  if [[ -f "${out}" && "${out}" -nt "${sam}" ]]; then
    continue
  fi
  echo "[stage-bams] sam->bam ${base}"
  samtools view -bS "${sam}" | samtools sort -o "${out}"
  samtools index "${out}"
done
ls -la parity/build/sam-indexed-bam
