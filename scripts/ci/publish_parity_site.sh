#!/usr/bin/env bash
# Update docs/parity-site/data from a workflow run and deploy to GitHub Pages (root).
#
# Env:
#   PARITY_SITE_SOURCE=nightly|genomewide
#   PARITY_SITE_JSON     (nightly: path to happy_summary.json)
#   PARITY_SITE_RUN_DIR  (genomewide: run directory)
#   GITHUB_SHA, GITHUB_SERVER_URL, GITHUB_REPOSITORY, GITHUB_RUN_ID
#   PARITY_SITE_DEPLOY=1  (default) use peaceiris-compatible publish dir staging
#   PARITY_SITE_COMMIT=1  commit data/ back to the current branch
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"

source="${PARITY_SITE_SOURCE:?set PARITY_SITE_SOURCE=nightly|genomewide}"
site_dir="${repo_root}/docs/parity-site"
run_url=""
if [[ -n "${GITHUB_SERVER_URL:-}" && -n "${GITHUB_REPOSITORY:-}" && -n "${GITHUB_RUN_ID:-}" ]]; then
  run_url="${GITHUB_SERVER_URL}/${GITHUB_REPOSITORY}/actions/runs/${GITHUB_RUN_ID}"
fi

args=(
  python3 "${repo_root}/scripts/parity/giab/update_public_dashboard.py"
  --source "${source}"
  --site-dir "${site_dir}"
  --commit-sha "${GITHUB_SHA:-unknown}"
  --workflow-run-url "${run_url}"
)

if [[ "${source}" == "nightly" ]]; then
  json="${PARITY_SITE_JSON:?set PARITY_SITE_JSON}"
  args+=(--json "${json}")
else
  run_dir="${PARITY_SITE_RUN_DIR:?set PARITY_SITE_RUN_DIR}"
  args+=(--run-dir "${run_dir}")
fi

"${args[@]}"

if [[ "${PARITY_SITE_COMMIT:-1}" == "1" ]]; then
  if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    git add docs/parity-site/data/history.json docs/parity-site/data/latest.json || true
    if ! git diff --cached --quiet; then
      git config user.name "${GIT_AUTHOR_NAME:-github-actions[bot]}" || true
      git config user.email "${GIT_AUTHOR_EMAIL:-41898282+github-actions[bot]@users.noreply.github.com}" || true
      git commit -m "docs(parity-site): update equivalence history (${GITHUB_SHA:-local})"
      if [[ "${PARITY_SITE_PUSH:-1}" == "1" ]]; then
        git push || echo "[parity-site] warning: git push failed (non-fatal)"
      fi
    else
      echo "[parity-site] no history.json changes to commit"
    fi
  fi
fi

echo "[parity-site] site ready at ${site_dir}"
