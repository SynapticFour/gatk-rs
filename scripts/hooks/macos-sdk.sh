#!/usr/bin/env bash
# Select a macOS SDK this machine's linker can read.
#
# macOS 27 Command Line Tools ship MacOSX27.sdk. libSystem.B.tbd in that
# SDK names the architecture arm64e.x1. ld-1267 rejects that name, so
# clang fails with "malformed file" / "unknown architecture" and every
# Rust build script fails before clippy or tests run. xcrun still picks
# that SDK because it matches the running OS, including when xcode-select
# points at an Xcode that only contains the 26 SDK.
#
# Decide from the .tbd text, not from a probe that can succeed while
# cargo's linker still opens the broken SDK. A caller-supplied SDKROOT
# is kept only when that SDK does not name arm64e.x1. On Linux this
# returns immediately. When a newer ld understands arm64e.x1, delete
# the marker from the check or the probe SDK will no longer contain it.

sdk_names_arm64e_x1() {
  local sdk="$1"
  local tbd
  [[ -n "${sdk}" && -d "${sdk}" ]] || return 1
  for tbd in "${sdk}/usr/lib/libSystem.B.tbd" "${sdk}/usr/lib/libSystem.tbd"; do
    if [[ -f "${tbd}" ]] && grep -q 'arm64e.x1' "${tbd}"; then
      return 0
    fi
  done
  return 1
}

use_linkable_macos_sdk() {
  if [[ "$(uname -s)" != "Darwin" ]]; then
    return 0
  fi

  local tmp active default_sdk resolved_default best best_ver
  local parent cand resolved ver seen tested newer
  active="${SDKROOT:-}"
  if [[ -z "${active}" ]]; then
    active="$(xcrun --show-sdk-path 2>/dev/null || true)"
  fi
  if ! sdk_names_arm64e_x1 "${active}"; then
    return 0
  fi

  tmp="$(mktemp -d)"
  default_sdk="${active}"
  resolved_default=""
  if [[ -d "${default_sdk}" ]]; then
    resolved_default="$(cd "${default_sdk}" && pwd -P)"
  fi

  best=""
  best_ver=""
  : >"${tmp}/parents"
  if [[ -n "${default_sdk}" && -d "${default_sdk}" ]]; then
    dirname "${default_sdk}" >>"${tmp}/parents"
  fi
  printf '%s\n' \
    /Library/Developer/CommandLineTools/SDKs \
    /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs \
    >>"${tmp}/parents"

  seen=" "
  tested=" "
  while IFS= read -r parent || [[ -n "${parent}" ]]; do
    [[ -d "${parent}" ]] || continue
    case "${seen}" in
      *" ${parent} "*) continue ;;
    esac
    seen="${seen}${parent} "
    for cand in "${parent}"/MacOSX*.sdk; do
      [[ -d "${cand}" ]] || continue
      resolved="$(cd "${cand}" && pwd -P)"
      if [[ -n "${resolved_default}" && "${resolved}" == "${resolved_default}" ]]; then
        continue
      fi
      if sdk_names_arm64e_x1 "${cand}"; then
        continue
      fi
      case "${tested}" in
        *" ${resolved} "*) continue ;;
      esac
      tested="${tested}${resolved} "
      if ! printf 'int main(void){return 0;}\n' | SDKROOT="${cand}" cc -x c - -o "${tmp}/b" >/dev/null 2>"${tmp}/b.err"; then
        continue
      fi
      ver="$(basename "${cand}" | sed -n 's/^MacOSX\([0-9][0-9.]*\).*/\1/p')"
      if [[ -z "${ver}" ]]; then
        ver="0"
      fi
      newer=0
      if [[ -z "${best}" ]]; then
        newer=1
      elif awk -v a="${ver}" -v b="${best_ver}" 'BEGIN {
        n = split(a, aa, ".")
        m = split(b, bb, ".")
        limit = n
        if (m > limit) limit = m
        for (i = 1; i <= limit; i++) {
          x = aa[i] + 0
          y = bb[i] + 0
          if (x > y) exit 0
          if (x < y) exit 1
        }
        exit 1
      }'; then
        newer=1
      fi
      if [[ "${newer}" -eq 1 ]]; then
        best="${cand}"
        best_ver="${ver}"
      fi
    done
  done <"${tmp}/parents"

  rm -rf "${tmp}"
  if [[ -z "${best}" ]]; then
    echo "pre-commit: ld cannot read the macOS 27 SDK (arm64e.x1) and no older SDK linked." >&2
    return 1
  fi

  export SDKROOT="${best}"
  echo "pre-commit: ld cannot read the macOS 27 SDK (arm64e.x1); using ${SDKROOT}" >&2
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
  set -euo pipefail
  use_linkable_macos_sdk
  printf 'SDKROOT=%s\n' "${SDKROOT-}"
fi
