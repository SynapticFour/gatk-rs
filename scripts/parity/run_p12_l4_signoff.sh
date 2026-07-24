#!/usr/bin/env bash
# P12 L4 production sign-off battery (FORMAT numeric parity on 66 Java-only sites).
#
# Usage:
#   export P12_REFERENCE="$PWD/parity/realworld/assets/hs37d5.simple.fa"
#   ./scripts/parity/run_p12_l4_signoff.sh
#
# Logs: parity/reports/p12_l4_signoff_<timestamp>.log
#
# L4.2 gate: p12_format_parity with P12_PHASE_E=1, P12_L4_JAVA_FORMAT unset (algorithmic).
# L4.1 harness: same test with P12_L4_JAVA_FORMAT=1 (fixture overlay at emit).
#
# Note: macOS BSD env(1) has no `env -u`; this script uses bash unset in subshells.
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"

if [[ -z "${P12_REFERENCE:-}" ]]; then
  export P12_REFERENCE="${repo_root}/parity/realworld/assets/hs37d5.simple.fa"
fi
if [[ ! -f "${P12_REFERENCE}" ]]; then
  echo "P12_REFERENCE not found: ${P12_REFERENCE}" >&2
  exit 1
fi

bam="${P12_BAM:-${repo_root}/parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam}"
if [[ ! -f "${bam}" ]]; then
  echo "P12 BAM not found: ${bam}" >&2
  echo "Stage with: ./scripts/parity/realworld/02_stage_na12878_20k_bam.sh" >&2
  exit 1
fi

stamp="$(date -u +%Y%m%dT%H%M%SZ)"
log="${repo_root}/parity/reports/p12_l4_signoff_${stamp}.log"
canonical="${repo_root}/parity/reports/p12_l4_signoff_canonical.log"
exec > >(tee -a "${log}") 2>&1

echo "=== P12 L4 sign-off ${stamp} ==="
echo "P12_REFERENCE=${P12_REFERENCE}"
echo "P12_BAM=${bam}"

# FORMAT/fixup modules are pub only under parity_harness; production binary stays feature-free.
PHC=(cargo test -p gatk-haplotypecaller --features parity_harness)

run_test() {
  local label="$1"
  shift
  echo ""
  echo "=== ${label} ==="
  "$@"
}

run_l4_harness() {
  (
    unset P12_L4_JAVA_FORMAT GATK_RS_P12_EVENT_REGISTRY
    export P12_PHASE_E=1
    export P12_REFERENCE="${P12_REFERENCE}"
    export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-${repo_root}/target}"
    "$@"
  )
}

run_l4_algorithmic() {
  (
    unset P12_L4_JAVA_FORMAT GATK_RS_P12_EVENT_REGISTRY
    export P12_PHASE_E=1
    export P12_REFERENCE="${P12_REFERENCE}"
    export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-${repo_root}/target}"
    "$@"
  )
}

run_l4_fixture_overlay() {
  (
    unset GATK_RS_P12_EVENT_REGISTRY
    export P12_PHASE_E=1
    export P12_L4_JAVA_FORMAT=1
    export P12_REFERENCE="${P12_REFERENCE}"
    export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-${repo_root}/target}"
    "$@"
  )
}

# --- fast (fixture contracts; no BAM timing) ---
run_test "L4.0 fixture contract (66 rows)" \
  "${PHC[@]}" --test p12_format_parity_test \
  p12_java_format_fixture_contract --release

run_test "L4.0 cluster FORMAT fixture (3 rows)" \
  "${PHC[@]}" --test p12_format_parity_test \
  p12_cluster_format_fixture --release

run_test "L4.1 emit fixup unit tests" \
  "${PHC[@]}" --test p12_java_format_fixup_test --release

# --- L3 regression (emit set must stay 66/66) ---
run_test "L3a emit gate (P12_PHASE_E=1)" \
  run_l4_harness "${PHC[@]}" --test p12_parity_gate_test \
  p12_parity_gate --release -- --ignored --nocapture

# --- L4 long gates (ignored; ~4–5 min each) ---
run_test "L4.2 algorithmic FORMAT (no fixture overlay)" \
  run_l4_algorithmic "${PHC[@]}" --test p12_format_parity_test \
  p12_format_parity --release -- --ignored --nocapture

run_test "L4.1 harness FORMAT (fixture overlay)" \
  run_l4_fixture_overlay "${PHC[@]}" --test p12_format_parity_test \
  p12_format_parity --release -- --ignored --nocapture

cp -f "${log}" "${canonical}"

echo ""
echo "=== P12 L4 sign-off PASS ==="
echo "log: ${log}"
echo "canonical: ${canonical}"
