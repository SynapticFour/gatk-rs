#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"

echo "[p5-haplotype-candidate-diff] Running frozen Java-export parity fixture checks"
cargo test -p gatk-haplotypecaller --test p5_haplotype_candidate_diff_test --locked
echo "[p5-haplotype-candidate-diff] passed"
