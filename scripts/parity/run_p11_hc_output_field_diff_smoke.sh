#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
report_dir="${repo_root}/parity/reports"
mkdir -p "${report_dir}"
cd "${repo_root}"

json_out="${report_dir}/p11_hc_output_field_diff_smoke.json"
java_out="${report_dir}/p11_hc_output_field_diff_smoke.java.vcf"
rust_out="${report_dir}/p11_hc_output_field_diff_smoke.rust.vcf"
java_bam="${report_dir}/p11_hc_output_field_diff_smoke.java_input.bam"
rust_log="${report_dir}/p11_hc_output_field_diff_smoke.rust.log"
target_dir="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target-parity}"
mkdir -p "${target_dir}"
reference="${repo_root}/parity/fixtures/p5_live_reference.fa"
synthetic_sam="${repo_root}/parity/fixtures/p11_java_positive.sam"
cache_bam="${repo_root}/parity/build/sam-indexed-bam/p11_java_positive.bam"

# Indexed BAM for both oracles — Rust HC requires a BAM/CRAM index for -L queries.
# Prefer CI/local samtools cache (no Docker), then samtools, then GATK SortSam.
ensure_indexed_bam() {
  if [[ -f "${cache_bam}" && ( -f "${cache_bam}.bai" || -f "${cache_bam%.bam}.bai" ) ]]; then
    java_bam="${cache_bam}"
    return 0
  fi
  if command -v samtools >/dev/null 2>&1; then
    mkdir -p "$(dirname "${java_bam}")"
    samtools view -bS "${synthetic_sam}" | samtools sort -o "${java_bam}"
    samtools index "${java_bam}"
    return 0
  fi
  docker run --rm --platform "${GATK_DOCKER_PLATFORM:-linux/amd64}" \
    -v "${repo_root}:${repo_root}" \
    -w "${repo_root}" \
    "${GATK_DOCKER_IMAGE:-us.gcr.io/broad-gatk/gatk:4.4.0.0}" \
    gatk SortSam \
    -I "${synthetic_sam}" \
    -O "${java_bam}" \
    -SO coordinate \
    --CREATE_INDEX true \
    --QUIET true >/dev/null 2>&1
}

ensure_indexed_bam

set +e
docker run --rm --platform "${GATK_DOCKER_PLATFORM:-linux/amd64}" \
  -v "${repo_root}:${repo_root}" \
  -w "${repo_root}" \
  "${GATK_DOCKER_IMAGE:-us.gcr.io/broad-gatk/gatk:4.4.0.0}" \
  gatk HaplotypeCaller \
  -R "${reference}" \
  -I "${java_bam}" \
  -O "${java_out}" \
  --standard-min-confidence-threshold-for-calling 0.0 \
  --verbosity ERROR >/dev/null 2>&1
java_exit=$?
set -e

set +e
CARGO_TARGET_DIR="${target_dir}" cargo run --quiet --bin gatk-rs -- \
  HaplotypeCaller \
  -R parity/fixtures/p5_live_reference.fa \
  -I "${java_bam}" \
  -O "${rust_out}" \
  -L chrLive:1-63 >"${rust_log}" 2>&1
rust_exit=$?
set -e

python3 - "${java_out}" "${rust_out}" "${json_out}" "${java_exit}" "${rust_exit}" "${repo_root}" <<'PY'
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(sys.argv[6]) / "scripts" / "parity"))
from p11_field_compare import (  # noqa: E402
    count_variants,
    first_variant_fields,
    smoke_status,
)

java = pathlib.Path(sys.argv[1])
rust = pathlib.Path(sys.argv[2])
out = pathlib.Path(sys.argv[3])
java_exit = int(sys.argv[4])
rust_exit = int(sys.argv[5])

java_variants = count_variants(java)
rust_variants = count_variants(rust)
java_first = first_variant_fields(java)
rust_first = first_variant_fields(rust)
status, notes, mismatches = smoke_status(
    java_exit=java_exit,
    java_variants=java_variants,
    rust_exit=rust_exit,
    rust_variants=rust_variants,
    java_first=java_first,
    rust_first=rust_first,
)

payload = {
    "label": "phase11-hc-output-field-diff-smoke",
    "status": status,
    "java_exit": java_exit,
    "rust_exit": rust_exit,
    "java_variant_record_count": java_variants,
    "rust_variant_record_count": rust_variants,
    "mismatches": mismatches,
    "java_first_variant": java_first,
    "rust_first_variant": rust_first,
    "notes": notes,
}
out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(
    f"[p11-field-diff] status={status} java_variants={java_variants} "
    f"rust_variants={rust_variants} rust_exit={rust_exit}"
)
# Soft-pass only when the Java oracle is unavailable (no Docker). Hard-fail otherwise.
if status in {"fail", "rust_fail", "divergent_activation", "pending_activation"}:
    raise SystemExit(1)
PY

# Keep logs/JSON for CI artifacts; drop bulky intermediates on success.
if [[ -f "${json_out}" ]]; then
  status="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["status"])' "${json_out}")"
  if [[ "${status}" == "pass" || "${status}" == "java_unavailable" ]]; then
    # Only remove report-dir BAM we created; never delete the shared sam-indexed cache.
    if [[ "${java_bam}" == "${report_dir}/"* ]]; then
      rm -f "${java_bam}" "${java_bam}.bai" "${java_bam%.bam}.bai"
    fi
    rm -f "${java_out}" "${rust_out}" "${java_out}.idx" "${rust_log}"
  fi
fi
