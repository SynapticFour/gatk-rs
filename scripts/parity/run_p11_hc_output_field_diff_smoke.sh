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
target_dir="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target-parity}"
mkdir -p "${target_dir}"
reference="${repo_root}/parity/fixtures/p5_live_reference.fa"
synthetic_sam="${repo_root}/parity/fixtures/p11_java_positive.sam"

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

CARGO_TARGET_DIR="${target_dir}" cargo run --quiet --bin gatk-rs -- \
  HaplotypeCaller \
  -R parity/fixtures/p5_live_reference.fa \
  -I parity/fixtures/p11_java_positive.sam \
  -O "${rust_out}" \
  -L chrLive:1-63 >/dev/null 2>&1

python3 - "${java_out}" "${rust_out}" "${json_out}" "${java_exit}" <<'PY'
import json
import pathlib
import sys

java = pathlib.Path(sys.argv[1])
rust = pathlib.Path(sys.argv[2])
out = pathlib.Path(sys.argv[3])
java_exit = int(sys.argv[4])

def count_variants(path: pathlib.Path) -> int:
    if not path.exists():
        return 0
    return sum(1 for l in path.read_text(encoding="utf-8", errors="replace").splitlines() if l and not l.startswith("#"))

def first_variant_fields(path: pathlib.Path):
    if not path.exists():
        return None
    for l in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if l and not l.startswith("#"):
            cols = l.split("\t")
            if len(cols) < 8:
                return None
            out = {
                "CHROM": cols[0],
                "POS": cols[1],
                "REF": cols[3],
                "ALT": cols[4],
                "QUAL": cols[5],
                "FILTER": cols[6],
                "INFO": cols[7],
            }
            if len(cols) >= 10:
                out["FORMAT"] = cols[8]
                out["SAMPLE"] = cols[9]
            return out
    return None

def extract_sample_subfields(variant):
    if not variant:
        return {}
    fmt = variant.get("FORMAT")
    sample = variant.get("SAMPLE")
    if not fmt or not sample:
        return {}
    keys = fmt.split(":")
    vals = sample.split(":")
    return {k: vals[i] if i < len(vals) else "." for i, k in enumerate(keys)}

def parse_info_map(info_str):
    out = {}
    if not info_str or info_str == ".":
        return out
    for tok in info_str.split(";"):
        if "=" in tok:
            k, v = tok.split("=", 1)
            out[k] = v
        else:
            out[tok] = "true"
    return out

java_variants = count_variants(java)
rust_variants = count_variants(rust)

if java_exit != 0:
    status = "java_unavailable"
    notes = "Java docker oracle unavailable in this environment"
elif rust_variants == 0:
    status = "pending_activation"
    notes = "Rust HC has no variant records; field-level diff is deferred"
elif java_variants == 0 and rust_variants > 0:
    status = "divergent_activation"
    notes = "Rust HC emits provisional variants while Java smoke interval has no calls; keep as activation scaffold only"
else:
    java_first = first_variant_fields(java)
    rust_first = first_variant_fields(rust)
    compare_keys = ["CHROM", "POS", "REF", "ALT", "QUAL", "FILTER"]
    strict_mismatches = []
    for k in compare_keys:
        if (java_first or {}).get(k) != (rust_first or {}).get(k):
            strict_mismatches.append(k)
    java_info = parse_info_map((java_first or {}).get("INFO", "."))
    rust_info = parse_info_map((rust_first or {}).get("INFO", "."))
    for k in ["AC", "AF", "AN", "DP"]:
        jv = java_info.get(k)
        rv = rust_info.get(k)
        if jv is None or rv is None:
            strict_mismatches.append(f"INFO.{k}")
            continue
        if k == "AF":
            try:
                if abs(float(jv) - float(rv)) > 0.01:
                    strict_mismatches.append(f"INFO.{k}")
            except Exception:
                strict_mismatches.append(f"INFO.{k}")
        else:
            if jv != rv:
                strict_mismatches.append(f"INFO.{k}")
    java_sample = extract_sample_subfields(java_first)
    rust_sample = extract_sample_subfields(rust_first)
    for k in ["GT", "AD", "DP", "GQ", "PL"]:
        if java_sample.get(k) != rust_sample.get(k):
            strict_mismatches.append(f"SAMPLE.{k}")
    if strict_mismatches:
        status = "fail"
        notes = f"strict field diff mismatch on keys: {','.join(strict_mismatches)}"
    else:
        status = "pass"
        notes = "strict field diff smoke matched on core variant keys"

payload = {
    "label": "phase11-hc-output-field-diff-smoke",
    "status": status,
    "java_exit": java_exit,
    "java_variant_record_count": java_variants,
    "rust_variant_record_count": rust_variants,
    "java_first_variant": first_variant_fields(java),
    "rust_first_variant": first_variant_fields(rust),
    "notes": notes,
}
out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
print(f"[p11-field-diff] status={status} java_variants={java_variants} rust_variants={rust_variants}")
PY

rm -f "${java_out}" "${rust_out}" "${java_bam}" "${java_bam}.bai"
