#!/usr/bin/env bash
# Phase B (B.1–B.4) — run all required L1 parity gates. Fails fast on first mismatch.
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

if [[ "${PARITY_SKIP_HC_FULL_PHASE_B:-0}" == "1" ]]; then
  echo "[hc-full-parity-phase-b] skipped (PARITY_SKIP_HC_FULL_PHASE_B=1)"
  exit 0
fi

echo "[hc-full-parity-phase-b] cargo test -p gatk-haplotypecaller --features parity_harness (skip unrelated Phase E assembler unit test)"
# Integration tests (e.g. p12_java_site_trace) need the pub P12 / discovery surface.
cargo test -p gatk-haplotypecaller --features parity_harness -- --skip p5_case1_assembler_emits_reference_haplotype --skip pairhmm_likelihood_vector_matches_frozen_java_dump_fixture

for gate in \
  run_hc_full_parity_b1_read_shards.sh \
  run_hc_full_parity_b2_locus.sh \
  run_hc_full_parity_b2_assembly_regions.sh \
  run_hc_full_parity_b3_apply_summary.sh \
  run_hc_full_parity_b4_walker_traversal.sh \
  run_hc_full_parity_b5_reads.sh \
  run_hc_full_parity_b5_ref.sh \
  run_hc_full_parity_b5_feature.sh \
  run_hc_full_parity_b5_trim.sh \
  run_hc_full_parity_b5_pileup_track.sh \
  run_hc_full_parity_b5_force_active.sh; do
  echo "[hc-full-parity-phase-b] ${gate}"
  "${repo_root}/scripts/parity/${gate}"
done

echo "[hc-full-parity-phase-b] Phase B L1 gates (B.1–B.4 + B.5) + unit tests: OK"
