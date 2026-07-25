#!/usr/bin/env bash

set -euo pipefail

MODE="${1:-report}"

# Scaffolding / experimental modules are compiled into gatk-core but are not the
# product surface for the line-coverage ratchet (would pin the gate near ~35%).
IGNORE_SCAFFOLD='(src/parallel/|src/benchmarking/)'

run_report() {
  cargo llvm-cov --package gatk-core --tests \
    --ignore-filename-regex "${IGNORE_SCAFFOLD}" \
    --summary-only
}

run_gate_minimum() {
  # Measured floor on gatk-core product modules after ignoring scaffolds (~44% lines
  # on 2026-07-25 main). Ratchet, do not invent coverage — raise only with real tests.
  cargo llvm-cov --package gatk-core --tests \
    --ignore-filename-regex "${IGNORE_SCAFFOLD}" \
    --fail-under-lines 40 \
    --summary-only
}

run_gate_priority() {
  # Phase C gate for prioritized modules. Keep this strict gate opt-in
  # until all unstable/generated areas are retired from active development.
  cargo llvm-cov \
    --package gatk-core \
    --tests \
    --fail-under-lines 90 \
    --include-path gatk-core/src/types.rs \
    --include-path gatk-core/src/utils.rs \
    --include-path gatk-core/src/io/fasta.rs \
    --include-path gatk-core/src/io/fastq.rs \
    --summary-only
}

case "${MODE}" in
  report)
    run_report
    ;;
  gate-minimum)
    run_gate_minimum
    ;;
  gate-priority)
    run_gate_priority
    ;;
  *)
    echo "Usage: $0 [report|gate-minimum|gate-priority]" >&2
    exit 1
    ;;
esac
