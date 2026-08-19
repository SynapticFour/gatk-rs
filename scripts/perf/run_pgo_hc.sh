#!/usr/bin/env bash
# Two-stage PGO for gatk-cli HaplotypeCaller (opt-in; NOT default release).
#
# Does not modify Cargo.toml. Requires llvm-profdata on PATH (rustup component
# llvm-tools-preview or system LLVM).
#
# Usage:
#   ./scripts/perf/run_pgo_hc.sh prepare   # instrumented build
#   ./scripts/perf/run_pgo_hc.sh train     # run representative HC workloads
#   ./scripts/perf/run_pgo_hc.sh optimize  # merge + PGO rebuild
#   ./scripts/perf/run_pgo_hc.sh compare   # print binary paths for fair A/B
#
# Env:
#   PGO_DIR           — working directory (default docs/perf/runs/pgo_<stamp>)
#   PGO_TRAIN_CMD     — shell command that runs HC (required for train)
#   PGO_SKIP_NATIVE=1 — do not add target-cpu=native (default: portable)
#
# Prove benefit with scripts/perf/run_fair_hc_comparison.sh against a non-PGO
# release binary before considering any default enablement.
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"
cd "${repo_root}"

cmd="${1:-}"
stamp="$(date -u +%Y%m%dT%H%M%SZ)"
PGO_DIR="${PGO_DIR:-${repo_root}/docs/perf/runs/pgo_${stamp}}"
mkdir -p "${PGO_DIR}"/{instrumented,merged,bin}
export CARGO_TARGET_DIR="${PGO_DIR}/target"

rustflags_base="${RUSTFLAGS:-}"
if [[ "${PGO_SKIP_NATIVE:-1}" != "1" ]]; then
  rustflags_base="${rustflags_base} -C target-cpu=native"
fi

need_profdata() {
  if ! command -v llvm-profdata >/dev/null 2>&1; then
    echo "llvm-profdata not found. Try: rustup component add llvm-tools-preview" >&2
    echo "Then: export PATH=\"\$(rustc --print sysroot)/lib/rustlib/\$(rustc -vV | awk '/host:/{print \$2}')/bin:\$PATH\"" >&2
    exit 1
  fi
}

case "${cmd}" in
  prepare)
    echo "[pgo] instrumented build → ${PGO_DIR}"
    RUSTFLAGS="${rustflags_base} -Cprofile-generate=${PGO_DIR}/instrumented" \
      cargo build -p gatk-cli --release --locked
    cp -f "${CARGO_TARGET_DIR}/release/gatk-rs" "${PGO_DIR}/bin/gatk-rs.instrumented" || \
      cp -f "${CARGO_TARGET_DIR}/release/gatk-cli" "${PGO_DIR}/bin/gatk-rs.instrumented"
    echo "[pgo] wrote ${PGO_DIR}/bin/gatk-rs.instrumented"
    ;;
  train)
    if [[ -z "${PGO_TRAIN_CMD:-}" ]]; then
      cat >&2 <<EOF
Set PGO_TRAIN_CMD to a representative HaplotypeCaller invocation, e.g.:

  export PGO_TRAIN_CMD='path/to/gatk-rs HaplotypeCaller -R … -I … -O … -L chr21:…'

Use the same windows as wall-losers / fair HC comparison when possible.
EOF
      exit 1
    fi
    echo "[pgo] training with PGO_TRAIN_CMD"
    # Prefer instrumented binary if present.
    if [[ -x "${PGO_DIR}/bin/gatk-rs.instrumented" ]]; then
      export PATH="${PGO_DIR}/bin:${PATH}"
    fi
    bash -c "${PGO_TRAIN_CMD}"
    echo "[pgo] raw profiles in ${PGO_DIR}/instrumented"
    ;;
  optimize)
    need_profdata
    llvm-profdata merge -o "${PGO_DIR}/merged/hc.profdata" "${PGO_DIR}/instrumented"
    echo "[pgo] merged ${PGO_DIR}/merged/hc.profdata"
    RUSTFLAGS="${rustflags_base} -Cprofile-use=${PGO_DIR}/merged/hc.profdata -Cllvm-args=-pgo-warn-mismatch" \
      cargo build -p gatk-cli --release --locked
    cp -f "${CARGO_TARGET_DIR}/release/gatk-rs" "${PGO_DIR}/bin/gatk-rs.pgo" 2>/dev/null || \
      cp -f "${CARGO_TARGET_DIR}/release/gatk-cli" "${PGO_DIR}/bin/gatk-rs.pgo"
    echo "[pgo] wrote ${PGO_DIR}/bin/gatk-rs.pgo"
    echo "[pgo] Compare with a clean release binary via fair HC wall before enabling by default."
    ;;
  compare)
    echo "Instrumented: ${PGO_DIR}/bin/gatk-rs.instrumented"
    echo "PGO optimized: ${PGO_DIR}/bin/gatk-rs.pgo"
    echo "Baseline: build with plain 'cargo build -p gatk-cli --release --locked' (separate target dir)."
    ;;
  *)
    echo "Usage: $0 {prepare|train|optimize|compare}" >&2
    exit 1
    ;;
esac
