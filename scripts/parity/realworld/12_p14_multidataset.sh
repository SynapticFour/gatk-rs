#!/usr/bin/env bash
# Optional — Phase 14 consolidated multi-dataset report (NA12878+GIAB + pending slots).
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "${repo_root}"
echo "=== 12_p14_multidataset ==="
./scripts/parity/run_p14_multidataset_equivalence.sh
echo "12_p14_multidataset: OK"
