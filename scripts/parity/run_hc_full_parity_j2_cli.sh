#!/usr/bin/env bash
# Phase J.2.1 — default CLI uses assembly-region-v1 (not provisional-output-v1).
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

tmp_vcf="${repo_root}/parity/reports/hc-full-parity-j2-cli-tmp.vcf"
mkdir -p "$(dirname "${tmp_vcf}")"

CARGO_TARGET_DIR="${CARGO_TARGET_DIR}" cargo run -q --bin gatk-rs -- \
  HaplotypeCaller \
  -R parity/fixtures/reference.fa \
  -I parity/fixtures/sample.bam \
  -O "${tmp_vcf}" \
  -L chr1:1-32

if ! grep -q 'assembly-region-v1' "${tmp_vcf}"; then
  echo "[hc-full-parity-j2-cli] expected GATK_RS_HC_PIPELINE=assembly-region-v1 in header" >&2
  grep 'GATK_RS_HC_PIPELINE' "${tmp_vcf}" >&2 || true
  exit 1
fi
if grep -q 'provisional-output-v1' "${tmp_vcf}"; then
  echo "[hc-full-parity-j2-cli] provisional-output-v1 must not be default pipeline" >&2
  exit 1
fi
echo "[hc-full-parity-j2-cli] ok assembly-region-v1 default pipeline"
rm -f "${tmp_vcf}"
