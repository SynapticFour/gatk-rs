#!/usr/bin/env bash
# Foundation layer — haplotypecaller library contracts (read model, traversal scaffolds).
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "${repo_root}"
echo "=== 10_test_haplotypecaller_lib ==="
cargo test -p gatk-haplotypecaller --lib --locked
echo "10_test_haplotypecaller_lib: OK"
