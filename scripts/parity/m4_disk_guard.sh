#!/usr/bin/env bash
# Lightweight disk/memory guard for M4 16GB runs.
# Usage:
#   source scripts/parity/m4_disk_guard.sh
#   m4_require_free_gb 8
#   m4_disk_watch 30 &  WATCH_PID=$!
#   ... work ...
#   kill $WATCH_PID 2>/dev/null || true
set -euo pipefail

m4_avail_gb() {
  # macOS: df -g reports GiB in column 4. GNU/Linux: df -Bg (or -BG) reports "123G".
  if df -g / >/dev/null 2>&1; then
    df -g / | awk 'NR==2{print int($4)}'
    return 0
  fi
  if df -Bg / >/dev/null 2>&1; then
    df -Bg / | awk 'NR==2{gsub(/[^0-9].*/, "", $4); print int($4)}'
    return 0
  fi
  # Fallback: 1K-blocks → GiB
  df -k / | awk 'NR==2{print int($4 / 1024 / 1024)}'
}

m4_require_free_gb() {
  local need="${1:-8}"
  local avail
  avail="$(m4_avail_gb)"
  if (( avail < need )); then
    echo "[m4-guard] FAIL: only ${avail}GB free (need ≥${need}GB)" >&2
    return 1
  fi
  echo "[m4-guard] ok: ${avail}GB free (need ≥${need}GB)"
}

m4_disk_watch() {
  local interval="${1:-60}"
  local min_gb="${2:-6}"
  while true; do
    local avail
    avail="$(m4_avail_gb)"
    local ts
    ts="$(date -u +%H:%M:%S)"
    echo "[m4-guard ${ts}] avail_GB=${avail}"
    if (( avail < min_gb )); then
      echo "[m4-guard] CRITICAL: avail_GB=${avail} < ${min_gb} — aborting watchers" >&2
      return 2
    fi
    sleep "$interval"
  done
}

# Prefer single workspace target; never spawn sandbox cargo-target dupes.
export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/target}"
export CARGO_BUILD_JOBS="${CARGO_BUILD_JOBS:-1}"
export CARGO_INCREMENTAL="${CARGO_INCREMENTAL:-0}"
# HC / rayon: keep RSS flat on 16GB machines (override explicitly if needed).
export RAYON_NUM_THREADS="${RAYON_NUM_THREADS:-2}"
# Release builds: skip fat LTO unless GATK_RS_RELEASE_LTO=1 (LTO OOMs / fills disk on M4).
if [[ "${GATK_RS_RELEASE_LTO:-0}" != "1" ]]; then
  export CARGO_PROFILE_RELEASE_LTO="${CARGO_PROFILE_RELEASE_LTO:-false}"
  export CARGO_PROFILE_RELEASE_CODEGEN_UNITS="${CARGO_PROFILE_RELEASE_CODEGEN_UNITS:-16}"
fi
