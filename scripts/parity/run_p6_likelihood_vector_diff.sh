#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"

echo "[p6-likelihood-vector-diff] running frozen likelihood-vector diff"
cargo test -p gatk-haplotypecaller --test p6_likelihood_vector_diff_test --locked
echo "[p6-likelihood-vector-diff] passed"
