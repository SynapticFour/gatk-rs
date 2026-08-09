#!/usr/bin/env bash
# Point git at the versioned hooks under .githooks/ (no per-clone copy into .git/hooks).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

if [[ ! -d .git ]]; then
  echo "error: not a git repository (run git init first)" >&2
  exit 1
fi

chmod +x .githooks/pre-commit .githooks/pre-push 2>/dev/null || true

git config core.hooksPath .githooks
echo "core.hooksPath=$(git config --get core.hooksPath)"
echo "Hooks installed:"
echo "  pre-commit — fmt, clippy (affected crates), ratchets, doc links, size, lib tests"
echo "  pre-push   — fmt + cargo clippy --workspace (matches CI lint gate)"
echo "Do not use --no-verify unless intentional; CI will still enforce the same checks."
