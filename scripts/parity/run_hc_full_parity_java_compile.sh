#!/usr/bin/env bash
# Compile HcFullParityGateDump against pinned GATK jar (Docker).
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=lib_pinned_gatk.sh
source "${script_dir}/lib_pinned_gatk.sh"

repo_root="${GATK_RS_REPO_ROOT}"
java_src="${repo_root}/scripts/parity/java/HcFullParityGateDump.java"
class_dir="${repo_root}/parity/build/hc-full-parity-java-classes"
docker_repo="${PARITY_HOST_REPO_ROOT:-${repo_root}}"

mkdir -p "${class_dir}"

if [[ "${PARITY_SKIP_JAVA_COMPILE:-0}" == "1" ]]; then
  echo "[hc-full-parity-java-compile] skipped (PARITY_SKIP_JAVA_COMPILE=1)"
  exit 0
fi

if ! command -v docker >/dev/null 2>&1; then
  echo "[hc-full-parity-java-compile] docker not available" >&2
  exit 127
fi

jar_name="${GATK_JAR_BASENAME:-gatk-package-4.4.0.0-local.jar}"
src_container="${docker_repo}/scripts/parity/java/HcFullParityGateDump.java"
out_container="${docker_repo}/parity/build/hc-full-parity-java-classes"

docker_args=(run --rm)
if [[ -n "${GATK_DOCKER_PLATFORM:-}" ]]; then
  docker_args+=(--platform "${GATK_DOCKER_PLATFORM}")
fi

echo "[hc-full-parity-java-compile] javac against ${GATK_DOCKER_IMAGE}"
docker "${docker_args[@]}" \
  -v "${docker_repo}:${docker_repo}" \
  -w "${docker_repo}" \
  "${GATK_DOCKER_IMAGE}" \
  bash -c "set -euo pipefail; mapfile -t src < <(find scripts/parity/java -name '*.java' | sort); javac -cp /gatk/${jar_name} -d '${out_container}' \"\${src[@]}\""

echo "[hc-full-parity-java-compile] classes -> ${class_dir}"
