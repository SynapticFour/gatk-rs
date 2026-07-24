#!/usr/bin/env bash
# Convenience wrapper: verify Phase-7 PL/GQ/AD/DP oracle + Phase-8 gVCF block oracle against frozen expected rows.

set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"

"${repo_root}/scripts/parity/run_p7_java_genotype_refresh.sh"
"${repo_root}/scripts/parity/run_p8_java_block_refresh.sh"

echo "[p78-java-refresh-bundle] both Java oracles matched frozen fixtures"
