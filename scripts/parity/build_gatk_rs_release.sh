#!/usr/bin/env bash
# Build gatk-rs release binary with visible progress (avoid silent cargo -q / IDE "Warming up..").
# Tuned for MacBook Air M4 16GB / limited free disk: -j1, no fat LTO by default.
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
# shellcheck source=m4_disk_guard.sh
source "${repo_root}/scripts/parity/m4_disk_guard.sh"
m4_require_free_gb "${GATK_RS_BUILD_MIN_FREE_GB:-8}" || exit 1
target_dir="${CARGO_TARGET_DIR:-${repo_root}/target}"
cd "${repo_root}"
jobs="${CARGO_BUILD_JOBS:-1}"
echo "[build-gatk-rs] target=${target_dir}/release/gatk-rs"
echo "[build-gatk-rs] jobs=${jobs} LTO=${CARGO_PROFILE_RELEASE_LTO:-workspace-default}"
echo "[build-gatk-rs] compiling gatk-haplotypecaller + gatk-cli (release)…"
CARGO_TARGET_DIR="${target_dir}" cargo build -p gatk-cli --release --bin gatk-rs -j "${jobs}"
echo "[build-gatk-rs] done: ${target_dir}/release/gatk-rs"
