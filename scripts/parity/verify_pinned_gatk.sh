#!/usr/bin/env bash
# P0 gate: local Java GATK build matches docs/GATK_PINNED.env.
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=lib_pinned_gatk.sh
source "${script_dir}/lib_pinned_gatk.sh"

repo_root="${GATK_RS_REPO_ROOT}"
pinned_env="${repo_root}/docs/GATK_PINNED.env"
sha_file="${repo_root}/GATK_PINNED_SHA"

if [[ "${PARITY_SKIP_GATK_PIN_VERIFY:-0}" == "1" ]]; then
  echo "[verify-pinned-gatk] skipped (PARITY_SKIP_GATK_PIN_VERIFY=1)"
  exit 0
fi

fail() {
  echo "[verify-pinned-gatk] FAIL: $*" >&2
  exit 1
}

[[ -f "${pinned_env}" ]] || fail "missing ${pinned_env}"
[[ -f "${sha_file}" ]] || fail "missing ${sha_file}"

file_sha="$(tr -d '[:space:]' <"${sha_file}")"
[[ "${file_sha}" == "${GATK_PINNED_SHA}" ]] || fail "GATK_PINNED_SHA file (${file_sha}) != env (${GATK_PINNED_SHA})"
[[ "${#file_sha}" -eq 40 ]] || fail "GATK_PINNED_SHA must be 40 hex chars, got ${#file_sha}"

version_ok() {
  local out="$1"
  grep -q "v${GATK_PINNED_REF}" <<<"${out}" || grep -q "${GATK_PINNED_REF}" <<<"${out}"
}

checked=0

if [[ -n "${GATK_JAR:-}" && -f "${GATK_JAR}" ]]; then
  out="$(java -jar "${GATK_JAR}" --version 2>&1 || true)"
  version_ok "${out}" || fail "GATK_JAR --version does not match ref ${GATK_PINNED_REF}: ${out}"
  echo "[verify-pinned-gatk] OK GATK_JAR=${GATK_JAR}"
  checked=1
elif command -v gatk >/dev/null 2>&1; then
  out="$(gatk --version 2>&1 || true)"
  version_ok "${out}" || fail "gatk on PATH --version does not match ref ${GATK_PINNED_REF}: ${out}"
  echo "[verify-pinned-gatk] OK gatk on PATH"
  checked=1
fi

if command -v docker >/dev/null 2>&1; then
  if ! docker image inspect "${GATK_DOCKER_IMAGE}" >/dev/null 2>&1; then
    echo "[verify-pinned-gatk] pulling ${GATK_DOCKER_IMAGE} …"
    if [[ -n "${GATK_DOCKER_PLATFORM:-}" ]]; then
      docker pull --platform "${GATK_DOCKER_PLATFORM}" "${GATK_DOCKER_IMAGE}"
    else
      docker pull "${GATK_DOCKER_IMAGE}"
    fi
  fi
  docker_args=(run --rm)
  if [[ -n "${GATK_DOCKER_PLATFORM:-}" ]]; then
    docker_args+=(--platform "${GATK_DOCKER_PLATFORM}")
  fi
  out="$(docker "${docker_args[@]}" "${GATK_DOCKER_IMAGE}" gatk --version 2>&1 || true)"
  version_ok "${out}" || fail "docker ${GATK_DOCKER_IMAGE} --version does not match ref ${GATK_PINNED_REF}: ${out}"
  echo "[verify-pinned-gatk] OK docker image ${GATK_DOCKER_IMAGE}"
  checked=1
fi

if [[ "${checked}" -eq 0 ]]; then
  echo "[verify-pinned-gatk] WARN: no GATK_JAR, gatk on PATH, or docker — pin files only (${GATK_PINNED_REF} @ ${GATK_PINNED_SHA_SHORT})"
  exit 0
fi

echo "[verify-pinned-gatk] pin OK: ${GATK_PINNED_REF} @ ${GATK_PINNED_SHA}"
