#!/usr/bin/env bash
# Phase 9 (113): Java GATK HaplotypeCaller smoke on parity fixtures (exit + parsable VCF only).

set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
report_dir="${repo_root}/parity/reports"
mkdir -p "${report_dir}"

gatk_image="${GATK_DOCKER_IMAGE:-us.gcr.io/broad-gatk/gatk:4.4.0.0}"
gatk_platform="${GATK_DOCKER_PLATFORM:-linux/amd64}"

out="${report_dir}/p9_java_hc_smoke.vcf"
stdout="${report_dir}/p9_java_hc_smoke.stdout.txt"

echo "[p9-java-hc-smoke] image=${gatk_image} platform=${gatk_platform}"

docker run --rm --platform "${gatk_platform}" \
  -v "${repo_root}:${repo_root}" \
  -w "${repo_root}" \
  "${gatk_image}" \
  gatk HaplotypeCaller \
  -R "${repo_root}/parity/fixtures/reference.fa" \
  -I "${repo_root}/parity/fixtures/sample.bam" \
  -O "${out}" \
  -L chr1:1-32 \
  --verbosity ERROR \
  >"${stdout}" 2>&1

if [[ ! -f "${out}" ]]; then
  echo "[p9-java-hc-smoke] missing output VCF at ${out}" >&2
  exit 1
fi

python3 - "${repo_root}" <<'PY'
from pathlib import Path
import sys

repo = Path(sys.argv[1])
vcf = repo / "parity" / "reports" / "p9_java_hc_smoke.vcf"
text = vcf.read_text(encoding="utf-8", errors="replace")
if "##fileformat=" not in text or "#CHROM" not in text:
    print("[p9-java-hc-smoke] java output does not look like a VCF", file=sys.stderr)
    raise SystemExit(1)
print("[p9-java-hc-smoke] java produced a parsable VCF header")
PY

rm -f "${out}" "${stdout}"

echo "[p9-java-hc-smoke] passed"
