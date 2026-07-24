#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"

echo "[p6-wave-c] step83/84 boundary-artifact contracts"
cargo test -p gatk-haplotypecaller --test p6_boundary_artifact_contract_test --locked

echo "[p6-wave-c] step86 fp policy contracts"
cargo test -p gatk-haplotypecaller --test p6_fp_policy_contract_test --locked

echo "[p6-wave-c] step85 pairhmm bench smoke"
cargo bench -p gatk-haplotypecaller --bench pairhmm --locked -- --quick

echo "[p6-wave-c] passed"
