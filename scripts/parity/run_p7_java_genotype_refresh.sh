#!/usr/bin/env bash
# Regenerate / verify parity/expected/p7_genotype_fields_smoke.java.tsv using the standalone Java oracle
# (same rounding rules as Rust emit_genotype_format_fields).
#
# Requires Docker + JDK inside the GATK image.
#
# Usage:
#   ./scripts/parity/run_p7_java_genotype_refresh.sh
#   WRITE_EXPECTED=1 ./scripts/parity/run_p7_java_genotype_refresh.sh

set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
expected="${repo_root}/parity/expected/p7_genotype_fields_smoke.java.tsv"

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
  bash -lc 'javac -encoding UTF-8 -d /tmp/p7parity /work/scripts/parity/java/P7GenotypeFieldsDump.java && java -cp /tmp/p7parity P7GenotypeFieldsDump /work/parity/fixtures/p7_genotype_fields_smoke.tsv' \
  >"${gen}"

if [[ "${WRITE_EXPECTED:-0}" == "1" ]]; then
  cp "${gen}" "${expected}"
  echo "[p7-java-refresh] wrote ${expected}"
else
  if ! diff -u "${expected}" "${gen}"; then
    echo "[p7-java-refresh] mismatch (set WRITE_EXPECTED=1 after review)" >&2
    exit 1
  fi
  echo "[p7-java-refresh] frozen expected matches Java oracle"
fi
