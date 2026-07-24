#!/usr/bin/env bash
# Foundation layer — gatk-core unit contracts (reference, intervals, cache).
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "${repo_root}"
echo "=== 08_test_gatk_core_lib ==="
cargo test -p gatk-core --lib --locked
echo "08_test_gatk_core_lib: OK"
