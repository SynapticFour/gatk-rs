#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "usage: $0 <stdout-file> <args...>" >&2
  exit 2
fi

stdout_file="$1"
shift

profile="${PARITY_RUST_PROFILE:-dev}"
target_dir="${PARITY_CARGO_TARGET_DIR:-$(pwd)/target-parity}"
mkdir -p "${target_dir}"

if [[ "${profile}" == "release" ]]; then
  CARGO_TARGET_DIR="${target_dir}" cargo run --quiet --release --bin gatk-rs -- "$@" >"${stdout_file}" 2>&1
else
  CARGO_TARGET_DIR="${target_dir}" cargo run --quiet --bin gatk-rs -- "$@" >"${stdout_file}" 2>&1
fi
