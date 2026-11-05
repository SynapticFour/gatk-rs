#!/usr/bin/env bash
# Append one equivalence record to gh-pages data/history.json.
#
# Live history lives on the gh-pages branch, not on main. This script fetches
# that file, appends one record, and pushes only the updated history blobs.
# It does not push the workflow checkout branch.
#
# GitHub Actions concurrency groups do not cross workflow files, so
# giab-genomewide, genomewide-validation, nightly-equivalence, and benchmark
# can write gh-pages together. The history push is a fast-forward of a commit
# whose parent is the tip just fetched. A rejected push refetches and merges
# again, so one publisher cannot drop another publisher's row.
#
# After a successful update, PARITY_SITE_OMIT_FROM_SITE_DEPLOY removes
# history.json and latest.json from docs/parity-site. The following peaceiris
# deploy copies that directory with keep_files, and therefore leaves the
# gh-pages copies of those two files in place. Do not put them in
# peaceiris exclude_assets: that option deletes the destination files.
#
# Fallback, and only this fallback: when the gh-pages branch does not exist,
# or it exists without data/history.json, the checkout copy is the base.
# A missing checkout copy starts an empty history. A failed fetch or a
# malformed history file aborts and does not write.
#
# Env:
#   PARITY_SITE_SOURCE=nightly|genomewide
#   PARITY_SITE_WORKFLOW   caller workflow id (required)
#   PARITY_SITE_JSON       nightly happy_summary.json
#   PARITY_SITE_RUN_DIR    genomewide run directory
#   PARITY_SITE_PUSH_PAGES=1
#       Push the history commit to gh-pages. Does not push main.
#   PARITY_SITE_OMIT_FROM_SITE_DEPLOY=1
#       Remove the checkout copies after the gh-pages update succeeds.
#   PARITY_SITE_COMMIT=0   do not commit history onto the checkout branch
#   PARITY_SITE_PUSH=0     do not push the checkout branch
#   PARITY_SITE_GIT_REMOTE default origin (a remote name or a local path)
#   PARITY_SITE_REPO_ROOT  default: repository containing this script
#   PARITY_SITE_PAGES_HISTORY_STATUS / PARITY_SITE_PAGES_HISTORY_FILE
#       Test override. Skips git fetch.
set -euo pipefail

if [[ -n "${PARITY_SITE_REPO_ROOT:-}" ]]; then
  repo_root="${PARITY_SITE_REPO_ROOT}"
else
  repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
fi
cd "${repo_root}"

source="${PARITY_SITE_SOURCE:?set PARITY_SITE_SOURCE=nightly|genomewide}"
workflow="${PARITY_SITE_WORKFLOW:?set PARITY_SITE_WORKFLOW to the generating workflow id}"
remote="${PARITY_SITE_GIT_REMOTE:-origin}"
site_dir="${repo_root}/docs/parity-site"
pages_ref="refs/parity-pages/tip"

run_url=""
if [[ -n "${GITHUB_SERVER_URL:-}" && -n "${GITHUB_REPOSITORY:-}" && -n "${GITHUB_RUN_ID:-}" ]]; then
  run_url="${GITHUB_SERVER_URL}/${GITHUB_REPOSITORY}/actions/runs/${GITHUB_RUN_ID}"
fi

pages_history_file="$(mktemp)"
fetch_err="$(mktemp)"
fallback_copy="$(mktemp)"
cleanup() {
  rm -f "${pages_history_file}" "${fetch_err}" "${fallback_copy}"
}
trap cleanup EXIT

if [[ -f "${site_dir}/data/history.json" ]]; then
  cp "${site_dir}/data/history.json" "${fallback_copy}"
fi

resolve_pages_status() {
  if [[ -n "${PARITY_SITE_PAGES_HISTORY_STATUS:-}" ]]; then
    case "${PARITY_SITE_PAGES_HISTORY_STATUS}" in
      present|absent|unavailable) ;;
      *)
        echo "[parity-site] bad PARITY_SITE_PAGES_HISTORY_STATUS" >&2
        exit 2
        ;;
    esac
    if [[ "${PARITY_SITE_PAGES_HISTORY_STATUS}" == "present" ]]; then
      if [[ -z "${PARITY_SITE_PAGES_HISTORY_FILE:-}" || ! -f "${PARITY_SITE_PAGES_HISTORY_FILE}" ]]; then
        echo "[parity-site] present status requires PARITY_SITE_PAGES_HISTORY_FILE" >&2
        exit 2
      fi
      cp "${PARITY_SITE_PAGES_HISTORY_FILE}" "${pages_history_file}"
    fi
    printf '%s\n' "${PARITY_SITE_PAGES_HISTORY_STATUS}"
    return 0
  fi

  local rc=0
  git ls-remote --exit-code "${remote}" refs/heads/gh-pages >/dev/null 2>"${fetch_err}" || rc=$?
  if [[ "${rc}" -eq 2 ]]; then
    echo "[parity-site] no gh-pages branch; checkout history is the base" >&2
    printf '%s\n' "absent"
    return 0
  fi
  if [[ "${rc}" -ne 0 ]]; then
    echo "[parity-site] gh-pages lookup failed; not replacing history" >&2
    cat "${fetch_err}" >&2 || true
    printf '%s\n' "unavailable"
    return 0
  fi
  if ! git fetch --depth=1 "${remote}" "+refs/heads/gh-pages:${pages_ref}" >"${fetch_err}" 2>&1; then
    echo "[parity-site] gh-pages fetch failed; not replacing history" >&2
    cat "${fetch_err}" >&2 || true
    printf '%s\n' "unavailable"
    return 0
  fi
  if git cat-file -e "${pages_ref}:data/history.json" 2>/dev/null; then
    git show "${pages_ref}:data/history.json" > "${pages_history_file}"
    printf '%s\n' "present"
    return 0
  fi
  echo "[parity-site] gh-pages has no data/history.json; checkout history is the base" >&2
  printf '%s\n' "absent"
}

omit_from_site_deploy() {
  if [[ "${PARITY_SITE_OMIT_FROM_SITE_DEPLOY:-0}" != "1" ]]; then
    return 0
  fi
  rm -f "${site_dir}/data/history.json" "${site_dir}/data/latest.json"
  echo "[parity-site] omitted history.json from the site publish directory so the deploy keeps the gh-pages copy"
}

history_commit_tree() {
  local tree="$1"
  local parent="${2:-}"
  if [[ -n "${parent}" ]]; then
    GIT_AUTHOR_NAME='github-actions[bot]' \
    GIT_AUTHOR_EMAIL='41898282+github-actions[bot]@users.noreply.github.com' \
    GIT_COMMITTER_NAME='github-actions[bot]' \
    GIT_COMMITTER_EMAIL='41898282+github-actions[bot]@users.noreply.github.com' \
    git commit-tree "${tree}" -p "${parent}" \
      -m "parity-site: record equivalence history (${GITHUB_SHA:-local})"
  else
    GIT_AUTHOR_NAME='github-actions[bot]' \
    GIT_AUTHOR_EMAIL='41898282+github-actions[bot]@users.noreply.github.com' \
    GIT_COMMITTER_NAME='github-actions[bot]' \
    GIT_COMMITTER_EMAIL='41898282+github-actions[bot]@users.noreply.github.com' \
    git commit-tree "${tree}" \
      -m "parity-site: record equivalence history (${GITHUB_SHA:-local})"
  fi
}

publish_history_commit() {
  local src_dir="$1"
  local index tree commit parent=""
  index="$(mktemp)"
  if git rev-parse --verify --quiet "${pages_ref}" >/dev/null; then
    parent="$(git rev-parse "${pages_ref}")" || return 1
    GIT_INDEX_FILE="${index}" git read-tree "${pages_ref}" || return 1
  else
    GIT_INDEX_FILE="${index}" git read-tree --empty || return 1
  fi
  mkdir -p "${src_dir}/data"
  GIT_INDEX_FILE="${index}" git --work-tree="${src_dir}" add -- data/history.json data/latest.json || return 1
  tree="$(GIT_INDEX_FILE="${index}" git write-tree)" || return 1
  rm -f "${index}"
  commit="$(history_commit_tree "${tree}" "${parent}")" || return 1
  git push "${remote}" "${commit}:refs/heads/gh-pages"
}

merge_once() {
  local status="$1"
  local dest="$2"
  local -a args
  args=(
    python3 "${repo_root}/scripts/parity/giab/update_public_dashboard.py"
    --source "${source}"
    --workflow "${workflow}"
    --site-dir "${dest}"
    --commit-sha "${GITHUB_SHA:-unknown}"
    --workflow-run-url "${run_url}"
    --pages-history-status "${status}"
  )
  if [[ "${status}" == "present" ]]; then
    args+=(--pages-history "${pages_history_file}")
  fi
  if [[ -s "${fallback_copy}" ]]; then
    args+=(--fallback-history "${fallback_copy}")
  fi
  if [[ "${source}" == "nightly" ]]; then
    args+=(--json "${PARITY_SITE_JSON:?set PARITY_SITE_JSON}")
  else
    args+=(--run-dir "${PARITY_SITE_RUN_DIR:?set PARITY_SITE_RUN_DIR}")
  fi
  "${args[@]}"
}

attempt=1
while [[ "${attempt}" -le 5 ]]; do
  status="$(resolve_pages_status)"
  if [[ "${status}" == "unavailable" ]]; then
    echo "[parity-site] refusing to publish equivalence history" >&2
    exit 1
  fi
  out_dir="$(mktemp -d)"
  merge_once "${status}" "${out_dir}"
  if [[ ! -f "${out_dir}/data/history.json" ]]; then
    echo "[parity-site] history record already present; left the stored record unchanged"
    omit_from_site_deploy
    exit 0
  fi
  if [[ "${PARITY_SITE_PUSH_PAGES:-0}" != "1" ]]; then
    mkdir -p "${site_dir}/data"
    cp "${out_dir}/data/history.json" "${site_dir}/data/history.json"
    cp "${out_dir}/data/latest.json" "${site_dir}/data/latest.json"
    echo "[parity-site] updated checkout history only (PARITY_SITE_PUSH_PAGES is not 1)"
    break
  fi
  if publish_history_commit "${out_dir}"; then
    echo "[parity-site] pushed equivalence history to gh-pages (no push to the checkout branch)"
    omit_from_site_deploy
    break
  fi
  echo "[parity-site] gh-pages advanced; retrying history merge (${attempt}/5)" >&2
  attempt=$((attempt + 1))
  if [[ "${attempt}" -gt 5 ]]; then
    echo "[parity-site] could not fast-forward gh-pages history without dropping records" >&2
    exit 1
  fi
done

# Defaults are off. Callers that still want a local commit must opt in.
# Pushing the checkout branch is not part of the Pages publish.
if [[ "${PARITY_SITE_COMMIT:-0}" == "1" ]]; then
  if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    git add docs/parity-site/data/history.json docs/parity-site/data/latest.json || true
    if ! git diff --cached --quiet; then
      git config user.name "${GIT_AUTHOR_NAME:-github-actions[bot]}" || true
      git config user.email "${GIT_AUTHOR_EMAIL:-41898282+github-actions[bot]@users.noreply.github.com}" || true
      git commit -m "docs(parity-site): update equivalence history (${GITHUB_SHA:-local})"
      if [[ "${PARITY_SITE_PUSH:-0}" == "1" ]]; then
        git push || echo "[parity-site] warning: git push of the checkout branch failed (non-fatal)"
      fi
    else
      echo "[parity-site] no checkout history.json changes to commit"
    fi
  fi
fi

echo "[parity-site] site ready at ${site_dir}"
