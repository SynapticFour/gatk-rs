#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"

echo "[p7-genotyping] running Phase-7 genotyping contract tests"
CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-1}" RAYON_NUM_THREADS="${RAYON_NUM_THREADS:-1}" \
  cargo test -p gatk-haplotypecaller --lib --locked genotyping::tests::
echo "[p7-genotyping] passed"
