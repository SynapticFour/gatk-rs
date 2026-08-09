#!/usr/bin/env bash
# Point git at the versioned hooks under .githooks/ (no per-clone copy into .git/hooks).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$ROOT"

if [[ ! -d .git ]]; then
  echo "error: not a git repository (run git init first)" >&2
  exit 1
fi

chmod +x .githooks/pre-commit .githooks/pre-push \
  scripts/dev/check_l2_g2_subset_live.sh 2>/dev/null || true

git config core.hooksPath .githooks
echo "core.hooksPath=$(git config --get core.hooksPath)"
echo "Hooks installed:"
echo "  pre-commit — fmt, clippy (affected crates), ratchets, md/path doc links,"
echo "               rustdoc broken intra-doc links (affected crates), excellence N-1…N-7,"
echo "               size, lib tests"
echo "  pre-push   — fmt, clippy --workspace, excellence N-1…N-7,"
echo "               rustdoc broken intra-doc links (--workspace),"
echo "               L2 g2-subset-live when assembly/k-best files changed"
echo "                 (skip with PARITY_HOOK_SKIP_G2=1; needs staged fixture BAMs)"
echo "Do not use --no-verify unless intentional; CI will still enforce the same checks."
