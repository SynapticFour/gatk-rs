#!/usr/bin/env bash
# Regenerate / verify parity/expected/p8_gvcf_blocks_smoke.java.tsv using the standalone Java oracle
# (same semantics as Rust build_gvcf_blocks_with_semantics).
#
# Requires Docker + JDK inside the GATK image. Intended for periodic refresh / audits — not called by CI unless invoked explicitly.
#
# Usage:
#   ./scripts/parity/run_p8_java_block_refresh.sh           # diff-only (exit 1 on mismatch)
#   WRITE_EXPECTED=1 ./scripts/parity/run_p8_java_block_refresh.sh   # overwrite frozen expected from Java oracle

set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
fixture="${repo_root}/parity/fixtures/p8_gvcf_blocks_smoke.tsv"
expected="${repo_root}/parity/expected/p8_gvcf_blocks_smoke.java.tsv"

gatk_image="${GATK_DOCKER_IMAGE:-us.gcr.io/broad-gatk/gatk:4.4.0.0}"
gatk_platform="${GATK_DOCKER_PLATFORM:-linux/amd64}"

gen="$(mktemp)"
cleanup() {
  rm -f "${gen}"
}
trap cleanup EXIT

docker run --rm --platform "${gatk_platform}" \
  -v "${repo_root}:/work" \
  -w /work "${gatk_image}" \
  bash -lc 'javac -encoding UTF-8 -d /tmp/p8parity /work/scripts/parity/java/P8GvcfBlockDump.java && java -cp /tmp/p8parity P8GvcfBlockDump /work/parity/fixtures/p8_gvcf_blocks_smoke.tsv' \
  >"${gen}"

if [[ "${WRITE_EXPECTED:-0}" == "1" ]]; then
  cp "${gen}" "${expected}"
  echo "[p8-java-refresh] wrote ${expected}"
else
  if ! diff -u "${expected}" "${gen}"; then
    echo "[p8-java-refresh] mismatch: ${expected} differs from Java oracle output (set WRITE_EXPECTED=1 after review)" >&2
    exit 1
  fi
  echo "[p8-java-refresh] frozen expected matches Java oracle"
fi
