#!/usr/bin/env bash
# Real-world playbook — step 01: toolchain sanity (no downloads).
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "${repo_root}"

echo "=== 01_check_environment (repo: ${repo_root}) ==="

fail=0
if ! command -v cargo >/dev/null 2>&1; then
  echo "MISSING: cargo (install Rust toolchain)" >&2
  fail=1
else
  echo "OK: cargo $(cargo --version)"
fi

if ! command -v docker >/dev/null 2>&1; then
  echo "MISSING: docker (required for Java GATK + samtools faidx in asset pipeline)" >&2
  fail=1
else
  echo "OK: docker $(docker --version | head -1)"
fi

img="${GATK_DOCKER_IMAGE:-us.gcr.io/broad-gatk/gatk:4.4.0.0}"
plat="${GATK_DOCKER_PLATFORM:-linux/amd64}"
echo "Checking GATK image (platform ${plat}): ${img}"
if ! docker image inspect "${img}" >/dev/null 2>&1; then
  echo "NOTE: image not present locally — first HC/asset step will pull (large)." >&2
else
  echo "OK: GATK image present locally"
fi

if command -v samtools >/dev/null 2>&1; then
  echo "OK: samtools $(samtools --version 2>&1 | head -1) (optional for some comparators)"
else
  echo "NOTE: samtools not on PATH — CI uses PARITY_REQUIRE_SAMTOOLS=1; install if you need local BAM parity tools."
fi

if [[ "${fail}" -ne 0 ]]; then
  echo "01_check_environment: FAILED (fix items above)" >&2
  exit 1
fi
echo "01_check_environment: OK"
exit 0
