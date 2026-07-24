#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"

echo "[p6-pairhmm] running step77-79 contracts"
cargo test -p gatk-haplotypecaller --test p6_pairhmm_contract_test --locked
