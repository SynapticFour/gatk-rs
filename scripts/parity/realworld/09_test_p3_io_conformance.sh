#!/usr/bin/env bash
# Foundation layer — Phase 3 IO roundtrip / indexed query integration tests.
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "${repo_root}"
echo "=== 09_test_p3_io_conformance ==="
cargo test -p gatk-core --test p3_io_conformance_tests --locked
echo "09_test_p3_io_conformance: OK"
