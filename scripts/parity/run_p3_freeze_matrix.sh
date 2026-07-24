#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
report_dir="${repo_root}/parity/reports"
mkdir -p "${report_dir}"

export LC_ALL=C
export TZ=UTC
export RAYON_NUM_THREADS="${RAYON_NUM_THREADS:-1}"
export PARITY_RANDOM_SEED="${PARITY_RANDOM_SEED:-1337}"
export PYTHONHASHSEED="${PYTHONHASHSEED:-${PARITY_RANDOM_SEED}}"
export SOURCE_DATE_EPOCH="${SOURCE_DATE_EPOCH:-1700000000}"

run_step() {
  local label="$1"
  shift
  echo "[p3-freeze] ${label}"
  "$@"
}

# Step-50 freeze matrix: pin both smoke profiles and full required P3 contract surface.
run_step "smoke-profile-smoke" env PARITY_SMOKE_PROFILE=smoke ./scripts/parity/run_parity_smoke.sh
run_step "smoke-profile-extended" env PARITY_SMOKE_PROFILE=extended ./scripts/parity/run_parity_smoke.sh

run_step "phase3-io-conformance" cargo test -p gatk-core --test p3_io_conformance_tests --locked
run_step "phase3-region-query-runtime-diff" ./scripts/parity/run_p3_region_query_diff.sh
run_step "phase3-region-records-runtime-diff" ./scripts/parity/run_p3_region_records_diff.sh
run_step "phase3-malformed-corpus-diff" ./scripts/parity/run_p3_malformed_corpus_diff.sh
run_step "phase3-indexed-edge-query-runtime-diff" ./scripts/parity/run_p3_indexed_edge_query_diff.sh
run_step "phase3-unmapped-supplementary-runtime-diff" ./scripts/parity/run_p3_unmapped_supplementary_diff.sh
run_step "phase3-truncation-corruption-runtime-diff" ./scripts/parity/run_p3_truncation_corruption_diff.sh
run_step "phase3-cram-roundtrip-contract" cargo test -p gatk-core --test p3_io_conformance_tests --locked cram_roundtrip_with_reference_preserves_optional_tags_contract
run_step "phase3-header-canonicalization-contract" cargo test -p gatk-core --test p3_io_conformance_tests --locked header_canonical_hd_sq_rg_stable_across_htslib_roundtrip_contract

echo "P3 freeze matrix completed successfully."
