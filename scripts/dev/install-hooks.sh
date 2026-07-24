#!/usr/bin/env bash
# Point git at the versioned hooks under .githooks/ (no per-clone copy into .git/hooks).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

if [[ ! -d .git ]]; then
  echo "error: not a git repository (run git init first)" >&2
  exit 1
fi

if [[ ! -x .githooks/pre-commit ]]; then
  chmod +x .githooks/pre-commit
fi

git config core.hooksPath .githooks
echo "core.hooksPath=$(git config --get core.hooksPath)"
echo "Pre-commit hook installed. Next commit will run fmt/clippy/ratchets/size checks."
