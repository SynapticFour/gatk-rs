#!/usr/bin/env bash
set -euo pipefail

_script_dir="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=lib_pinned_gatk.sh
source "${_script_dir}/lib_pinned_gatk.sh"

if [[ $# -lt 2 ]]; then
  echo "usage: $0 <stdout-file> <args...>" >&2
  exit 2
fi

stdout_file="$1"
shift

repo_root="${GATK_RS_REPO_ROOT}"
docker_repo_root="${PARITY_HOST_REPO_ROOT:-${repo_root}}"

translated_args=()
for arg in "$@"; do
  if [[ "${arg}" == "${repo_root}"* ]]; then
    translated_args+=( "${docker_repo_root}${arg#${repo_root}}" )
  else
    translated_args+=( "${arg}" )
  fi
done

if [[ -n "${GATK_JAR:-}" ]]; then
  java -jar "${GATK_JAR}" "$@" >"${stdout_file}" 2>&1
elif command -v gatk >/dev/null 2>&1; then
  gatk "$@" >"${stdout_file}" 2>&1
elif [[ -n "${GATK_DOCKER_IMAGE:-}" ]]; then
  if [[ -n "${GATK_DOCKER_PLATFORM:-}" ]]; then
    docker run --rm \
      --platform "${GATK_DOCKER_PLATFORM}" \
      -v "${docker_repo_root}:${docker_repo_root}" \
      -w "${docker_repo_root}" \
      "${GATK_DOCKER_IMAGE}" \
      gatk "${translated_args[@]}" >"${stdout_file}" 2>&1
  else
    docker run --rm \
      -v "${docker_repo_root}:${docker_repo_root}" \
      -w "${docker_repo_root}" \
      "${GATK_DOCKER_IMAGE}" \
      gatk "${translated_args[@]}" >"${stdout_file}" 2>&1
  fi
else
  echo "Missing Java GATK: set GATK_JAR, install gatk on PATH, or set GATK_DOCKER_IMAGE" >"${stdout_file}"
  exit 127
fi
