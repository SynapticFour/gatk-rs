#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"

echo "[p5-assembly-stability] repeated-run determinism checks"
RAYON_NUM_THREADS=1 cargo test -p gatk-haplotypecaller --test p5_assembly_regression_test --locked outputs_are_stable_across_repeated_runs_and_input_order
RAYON_NUM_THREADS=4 cargo test -p gatk-haplotypecaller --test p5_assembly_regression_test --locked outputs_are_stable_across_repeated_runs_and_input_order
echo "[p5-assembly-stability] passed"
