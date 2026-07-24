# shellcheck shell=bash
# Load pinned GATK upstream defaults from docs/GATK_PINNED.env.
# Safe to source multiple times. Does not override variables already set in the environment.
#
# Usage:
#   source "$(dirname "$0")/lib_pinned_gatk.sh"

if [[ -n "${_GATK_RS_PINNED_GATK_LOADED:-}" ]]; then
  return 0 2>/dev/null || exit 0
fi
_GATK_RS_PINNED_GATK_LOADED=1

_lib_dir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
_repo_root="$(cd "${_lib_dir}/../.." && pwd)"
_pinned_env="${_repo_root}/docs/GATK_PINNED.env"
_pinned_sha_file="${_repo_root}/GATK_PINNED_SHA"

if [[ ! -f "${_pinned_env}" ]]; then
  echo "lib_pinned_gatk: missing ${_pinned_env}" >&2
  return 1 2>/dev/null || exit 1
fi

# Export only keys not already set (respect user overrides).
while IFS= read -r line || [[ -n "${line}" ]]; do
  [[ -z "${line}" || "${line}" =~ ^[[:space:]]*# ]] && continue
  if [[ "${line}" =~ ^([A-Za-z_][A-Za-z0-9_]*)=(.*)$ ]]; then
    key="${BASH_REMATCH[1]}"
    val="${BASH_REMATCH[2]}"
    if [[ -z "${!key:-}" ]]; then
      export "${key}=${val}"
    fi
  fi
done <"${_pinned_env}"

# Root one-liner must match env file.
if [[ -f "${_pinned_sha_file}" ]]; then
  _file_sha="$(tr -d '[:space:]' <"${_pinned_sha_file}")"
  if [[ -n "${GATK_PINNED_SHA:-}" && "${_file_sha}" != "${GATK_PINNED_SHA}" ]]; then
    echo "lib_pinned_gatk: GATK_PINNED_SHA mismatch between env and ${_pinned_sha_file}" >&2
    return 1 2>/dev/null || exit 1
  fi
  export GATK_PINNED_SHA="${GATK_PINNED_SHA:-${_file_sha}}"
fi

export GATK_DOCKER_IMAGE="${GATK_DOCKER_IMAGE:-us.gcr.io/broad-gatk/gatk:4.4.0.0}"
export GATK_DOCKER_PLATFORM="${GATK_DOCKER_PLATFORM:-linux/amd64}"

gatk_github_tree_url() {
  local relpath="${1:-}"
  relpath="${relpath#/}"
  echo "https://github.com/broadinstitute/gatk/blob/${GATK_PINNED_SHA}/${relpath}"
}

export GATK_RS_REPO_ROOT="${_repo_root}"
