#!/usr/bin/env bash
# Phase 9 (113): verify Rust HaplotypeCaller scaffold VCF matches frozen golden (deterministic header).

set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"

out="$(mktemp)"
cleanup() {
  rm -f "${out}"
}
trap cleanup EXIT

CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-1}" RAYON_NUM_THREADS="${RAYON_NUM_THREADS:-1}" \
  GATK_RS_HC_SCAFFOLD_OUTPUT=1 \
  cargo run -q --bin gatk-rs -- \
  HaplotypeCaller \
  -R parity/fixtures/reference.fa \
  -I parity/fixtures/sample.bam \
  -O "${out}" \
  -L chr1:1-32

if ! diff -u parity/expected/p9_hc_scaffold_golden.vcf "${out}"; then
  echo "[p9-hc-scaffold] output differs from golden (refresh golden only after intentional contract change)" >&2
  exit 1
fi
echo "[p9-hc-scaffold] rust scaffold VCF matches golden"
