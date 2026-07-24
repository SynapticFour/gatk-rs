#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"

echo "[p8-gvcf] running Phase-8 gVCF contract tests"

# Step 101-103 gVCF block construction + END boundary semantics
CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-1}" RAYON_NUM_THREADS="${RAYON_NUM_THREADS:-1}" \
  cargo test -p gatk-haplotypecaller --lib --locked genotyping::tests::gvcf_

# Step 105 mode-switch semantics
CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-1}" RAYON_NUM_THREADS="${RAYON_NUM_THREADS:-1}" \
  cargo test -p gatk-haplotypecaller --lib --locked genotyping::tests::emit_mode_

# Step 106 no-variation region behavior
CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-1}" RAYON_NUM_THREADS="${RAYON_NUM_THREADS:-1}" \
  cargo test -p gatk-haplotypecaller --lib --locked genotyping::tests::no_variation_region_

# Step 107 joint-compat sanity checks
CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-1}" RAYON_NUM_THREADS="${RAYON_NUM_THREADS:-1}" \
  cargo test -p gatk-haplotypecaller --lib --locked genotyping::tests::joint_compat_

echo "[p8-gvcf] passed"
