#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"

export CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-1}"
export RAYON_NUM_THREADS="${RAYON_NUM_THREADS:-1}"

echo "[p9-cli] gatk-cli integration tests"
cargo test -p gatk-cli --test haplotype_caller_cli_integration --locked

echo "[p9-cli] haplotypecaller scaffold warmup"
cargo test -p gatk-haplotypecaller --test run_warmup_test --locked

echo "[p9-cli] gatk_cli_exit_code unit tests"
cargo test -p gatk-common --lib --locked gatk_cli_exit_code_mappings

echo "[p9-cli] passed"
