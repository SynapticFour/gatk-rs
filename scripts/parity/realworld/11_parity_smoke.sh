#!/usr/bin/env bash
# Foundation layer — Java vs Rust differential smoke (Docker + fixtures).
# Default to the same Docker image as other parity scripts so Java is not missing (exit 127).
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "${repo_root}"
export GATK_DOCKER_IMAGE="${GATK_DOCKER_IMAGE:-us.gcr.io/broad-gatk/gatk:4.4.0.0}"
export GATK_DOCKER_PLATFORM="${GATK_DOCKER_PLATFORM:-linux/amd64}"
echo "=== 11_parity_smoke (PARITY_SMOKE_PROFILE=${PARITY_SMOKE_PROFILE:-smoke}, GATK_DOCKER_IMAGE=${GATK_DOCKER_IMAGE}) ==="
./scripts/parity/run_parity_smoke.sh
echo "11_parity_smoke: OK"
