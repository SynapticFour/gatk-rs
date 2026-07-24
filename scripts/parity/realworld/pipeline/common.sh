#!/usr/bin/env bash
# shellcheck source=common.sh
# Shared env for Real-World pipeline steps. Source from other scripts in this directory.
set -euo pipefail

_PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
# pipeline/ -> realworld/ -> parity/ -> scripts/ -> repo
export REPO_ROOT="$(cd "${_PIPELINE_DIR}/../../../.." && pwd)"

export RW_REF="${RW_REF:-${REPO_ROOT}/parity/realworld/assets/hs37d5.simple.fa}"
export RW_BAM="${RW_BAM:-${REPO_ROOT}/parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam}"
export RW_INTERVAL="${RW_INTERVAL:-20:413419-463418}"

export OUT_DIR="${OUT_DIR:-${REPO_ROOT}/parity/reports/realworld_pipeline_run}"
mkdir -p "${OUT_DIR}"

# shellcheck source=../../lib_pinned_gatk.sh
source "${REPO_ROOT}/scripts/parity/lib_pinned_gatk.sh"
export PARITY_RUST_PROFILE="${PARITY_RUST_PROFILE:-release}"
export PARITY_CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${REPO_ROOT}/target-parity}"

export RUN_JAVA="${REPO_ROOT}/scripts/parity/run_java_gatk.sh"
export RUN_RUST="${REPO_ROOT}/scripts/parity/run_rust_gatk.sh"
