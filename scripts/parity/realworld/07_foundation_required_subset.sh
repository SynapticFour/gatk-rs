#!/usr/bin/env bash
# Real-world playbook — step 07: full required foundation gate (long). Run after real-world slice
# when you want "everything else still green" before changing fixtures.
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "${repo_root}"

echo "=== 07_foundation_required_subset (run_foundation_gate.sh) ==="
./scripts/parity/run_foundation_gate.sh
echo "07: all required checks from parity/checks.json passed"
