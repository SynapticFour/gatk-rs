#!/usr/bin/env bash
# Free GitHub-hosted runner disk so rust-lld does not bus-error (signal 7) on huge
# --all-features / parity_harness test binaries. No-op-ish on macOS (brew paths differ).
set -euo pipefail

df -h || true
if [[ "$(uname -s)" == "Linux" ]]; then
  sudo rm -rf \
    /usr/share/dotnet \
    /usr/local/lib/android \
    /opt/ghc \
    /usr/local/share/powershell \
    /usr/share/swift \
    /opt/hostedtoolcache/CodeQL \
    /usr/local/.ghcup \
    || true
  sudo apt-get clean || true
fi
df -h || true
