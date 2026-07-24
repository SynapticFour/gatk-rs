#!/usr/bin/env bash
# P12 L3 production sign-off battery (emit parity: 66/66 Java sites, rust_only=0).
#
# Usage:
#   export P12_REFERENCE="$PWD/parity/realworld/assets/hs37d5.simple.fa"
#   ./scripts/parity/run_p12_l3_signoff.sh
#
# Logs: parity/reports/p12_l3_signoff_<timestamp>.log
#
# Note: macOS BSD env(1) has no `env -u`; this script uses bash unset in subshells.
# All HC integration tests need `--features parity_harness` (P12 probe modules are
# pub only under that feature; production binary remains feature-free).
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

stamp="$(date -u +%Y%m%dT%H%M%SZ)"
log="${repo_root}/parity/reports/p12_l3_signoff_${stamp}.log"
exec > >(tee -a "${log}") 2>&1

echo "=== P12 L3 sign-off ${stamp} ==="
echo "P12_REFERENCE=${P12_REFERENCE}"

PHC=(cargo test -p gatk-haplotypecaller --features parity_harness)

run_test() {
  local label="$1"
  shift
  echo ""
  echo "=== ${label} ==="
  "$@"
}

# Graph-only production default (bridges off, no P12_PHASE_E whitelist).
run_graph_only() {
  (
    unset P12_PHASE_E GATK_RS_ASM8_ONLY GATK_RS_P12_ENSURE_BRIDGES
    export P12_REFERENCE="${P12_REFERENCE}"
    "$@"
  )
}

# L3a harness emit whitelist.
run_l3a_harness() {
  (
    unset GATK_RS_ASM8_ONLY GATK_RS_P12_ENSURE_BRIDGES
    export P12_PHASE_E=1
    export P12_REFERENCE="${P12_REFERENCE}"
    "$@"
  )
}

# Explicit ASM-8 env smoke (redundant with graph-only default).
run_asm8_explicit() {
  (
    unset P12_PHASE_E GATK_RS_P12_ENSURE_BRIDGES
    export GATK_RS_ASM8_ONLY=1
    export P12_REFERENCE="${P12_REFERENCE}"
    "$@"
  )
}

# --- fast regression (unit + probes) ---
run_test "harness cfg audit" \
  python3 scripts/parity/p12_site_id_audit.py --check-harness

run_test "unit: bridge flags" \
  "${PHC[@]}" --test p12_asm8_only_gate_test \
  p12_asm8_only_bridge_flags --release

run_test "unit: strict defaults" \
  "${PHC[@]}" --test strict_event_map_test --release

run_test "unit: call_region phase b" \
  "${PHC[@]}" --test call_region_phase_b_test --release

run_test "probe: 92309492 haplotype-pair rust-only" \
  "${PHC[@]}" --test p12_site_92309492_probe_test --release

run_test "probe: mid-A gap 92316315 (graph-only)" \
  run_graph_only "${PHC[@]}" \
  --test p12_region_923162_mid_a_genotyping_test \
  p12_region_923162_mid_a_genotyping --release -- --nocapture

run_test "probe: mid-B 92318210/92318227 (graph-only)" \
  run_graph_only "${PHC[@]}" \
  --test p12_region_923181_genotyping_probe_test --release -- --nocapture

# --- long gates (ignored) ---
run_test "L3a harness (P12_PHASE_E=1)" \
  run_l3a_harness "${PHC[@]}" --test p12_parity_gate_test \
  p12_parity_gate --release -- --ignored --nocapture

run_test "L3b production (default graph-only)" \
  run_graph_only "${PHC[@]}" \
  --test p12_production_parity_gate_test \
  p12_production_parity_gate --release -- --ignored --nocapture

run_test "ASM-8 production gate" \
  run_asm8_explicit "${PHC[@]}" \
  --test p12_asm8_production_gate_test \
  p12_asm8_production_parity_gate --release -- --ignored --nocapture

run_test "site trace (66× per-site)" \
  run_graph_only "${PHC[@]}" \
  --test p12_java_site_trace_test \
  p12_java_site_trace --release -- --ignored --nocapture

echo ""
echo "=== P12 L3 sign-off PASS ==="
echo "log: ${log}"
