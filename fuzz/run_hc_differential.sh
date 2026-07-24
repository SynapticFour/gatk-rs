#!/usr/bin/env bash
# HC differential fuzzer entrypoint (M4-safe defaults).
#
# Two modes:
#   1) Campaign (default): synthetic BAM → Java GATK4 vs gatk-rs → shrink → fixture
#      Optional: --open-github-issue (requires authenticated `gh`)
#   2) LibFuzzer generative smoke: cargo +nightly fuzz run hc_differential
#
# Usage:
#   ./fuzz/run_hc_differential.sh
#   ./fuzz/run_hc_differential.sh --iterations 4 --open-github-issue
#   ./fuzz/run_hc_differential.sh --libfuzzer
#   ./fuzz/run_hc_differential.sh --replay-fixture gatk-haplotypecaller/tests/fixtures/regressions/<id>

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"

# shellcheck disable=SC1091
if [[ -f scripts/parity/m4_disk_guard.sh ]]; then
  # shellcheck source=/dev/null
  source scripts/parity/m4_disk_guard.sh || true
fi

export CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-1}"
export RAYON_NUM_THREADS="${RAYON_NUM_THREADS:-2}"

RUST_BIN="${RUST_BIN:-$ROOT/target/release/gatk-rs}"
JAVA_JAR="${JAVA_GATK_JAR:-}"
JAVA_BIN="${JAVA_GATK_BIN:-}"
ITERATIONS="${ITERATIONS:-8}"
OUT="${DIFF_FUZZ_OUT:-$ROOT/target/diff-fuzz}"
FIXTURE_ROOT="${FIXTURE_ROOT:-$ROOT/gatk-haplotypecaller/tests/fixtures/regressions}"
FORMAT_AD_TOL="${FORMAT_AD_TOL:-0}"
OPEN_ISSUE=0
LIBFUZZER=0
REPLAY=""
EXTRA=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --libfuzzer) LIBFUZZER=1; shift ;;
    --open-github-issue) OPEN_ISSUE=1; shift ;;
    --iterations) ITERATIONS="$2"; shift 2 ;;
    --rust-binary) RUST_BIN="$2"; shift 2 ;;
    --java-gatk-jar) JAVA_JAR="$2"; shift 2 ;;
    --java-gatk-bin) JAVA_BIN="$2"; shift 2 ;;
    --out) OUT="$2"; shift 2 ;;
    --fixture-root) FIXTURE_ROOT="$2"; shift 2 ;;
    --format-ad-tol) FORMAT_AD_TOL="$2"; shift 2 ;;
    --replay-fixture) REPLAY="$2"; shift 2 ;;
    *) EXTRA+=("$1"); shift ;;
  esac
done

if [[ "$LIBFUZZER" -eq 1 ]]; then
  echo "[diff-fuzz] libFuzzer generative target (scenario decode only)"
  if command -v cargo-fuzz >/dev/null 2>&1 || cargo fuzz --help >/dev/null 2>&1; then
    exec cargo +nightly fuzz run hc_differential -- "${EXTRA[@]+"${EXTRA[@]}"}"
  fi
  echo "[diff-fuzz] cargo-fuzz missing; building lean check target instead"
  exec cargo check --manifest-path "$ROOT/fuzz/Cargo.toml" --no-default-features --bin hc_differential
fi

if [[ ! -x "$RUST_BIN" && ! -f "$RUST_BIN" ]]; then
  echo "[diff-fuzz] building gatk-rs release (jobs=1)…"
  cargo build --release -p gatk-cli --bin gatk-rs
  RUST_BIN="$ROOT/target/release/gatk-rs"
fi

ARGS=(
  differential-fuzz
  --iterations "$ITERATIONS"
  --out "$OUT"
  --rust-binary "$RUST_BIN"
  --fixture-root "$FIXTURE_ROOT"
  --format-ad-tol "$FORMAT_AD_TOL"
  --shrink-steps "${SHRINK_STEPS:-24}"
  --min-free-gb "${MIN_FREE_GB:-8}"
)

if [[ -n "$JAVA_JAR" ]]; then
  ARGS+=(--java-gatk-jar "$JAVA_JAR")
fi
if [[ -n "$JAVA_BIN" ]]; then
  ARGS+=(--java-gatk-bin "$JAVA_BIN")
fi
if [[ "$OPEN_ISSUE" -eq 1 ]]; then
  ARGS+=(--open-github-issue)
fi
if [[ -n "$REPLAY" ]]; then
  ARGS+=(--replay-fixture "$REPLAY")
fi
ARGS+=("${EXTRA[@]+"${EXTRA[@]}"}")

echo "[diff-fuzz] cargo run -p gatk-rs-equiv -- ${ARGS[*]}"
exec cargo run -p gatk-rs-equiv -- "${ARGS[@]}"
