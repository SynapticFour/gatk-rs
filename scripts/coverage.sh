#!/usr/bin/env bash

set -euo pipefail

MODE="${1:-report}"

run_report() {
  cargo llvm-cov --package gatk-core --tests --summary-only
}

run_gate_minimum() {
  cargo llvm-cov --package gatk-core --tests --fail-under-lines 70 --summary-only
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
