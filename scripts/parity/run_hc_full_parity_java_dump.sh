#!/usr/bin/env bash
# Run HcFullParityGateDump subcommand (stdout TSV). Requires run_hc_full_parity_java_compile.sh.
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "usage: $0 <subcommand> [args...]" >&2
  exit 2
fi

script_dir="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=lib_pinned_gatk.sh
source "${script_dir}/lib_pinned_gatk.sh"

repo_root="${GATK_RS_REPO_ROOT}"
class_dir="${repo_root}/parity/build/hc-full-parity-java-classes"
docker_repo="${PARITY_HOST_REPO_ROOT:-${repo_root}}"
jar_name="${GATK_JAR_BASENAME:-gatk-package-4.4.0.0-local.jar}"

if [[ ! -d "${class_dir}" ]] || [[ ! -f "${class_dir}/HcFullParityGateDump.class" ]]; then
  "${script_dir}/run_hc_full_parity_java_compile.sh"
fi

if command -v docker >/dev/null 2>&1 && [[ -n "${GATK_DOCKER_IMAGE:-}" ]]; then
  docker_args=(run --rm)
  if [[ -n "${GATK_DOCKER_PLATFORM:-}" ]]; then
    docker_args+=(--platform "${GATK_DOCKER_PLATFORM}")
  fi
  docker "${docker_args[@]}" \
    -v "${docker_repo}:${docker_repo}" \
    -w "${docker_repo}" \
    "${GATK_DOCKER_IMAGE}" \
    java -cp "/gatk/${jar_name}:${docker_repo}/parity/build/hc-full-parity-java-classes" \
    HcFullParityGateDump "$@"
  exit $?
fi

if [[ -n "${GATK_JAR:-}" && -f "${GATK_JAR}" ]]; then
  java -cp "${GATK_JAR}:${class_dir}" HcFullParityGateDump "$@"
  exit $?
fi

echo "No Java runtime: need docker (${GATK_DOCKER_IMAGE}) or GATK_JAR" >&2
exit 127
