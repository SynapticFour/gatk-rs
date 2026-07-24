#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
report_dir="${repo_root}/parity/reports"
mkdir -p "${report_dir}"
cd "${repo_root}"

tmp_vcf="${report_dir}/p11_hc_activation_contracts.vcf"
json_out="${report_dir}/p11_hc_output_activation_contracts.json"
target_dir="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target-parity}"
mkdir -p "${target_dir}"

CARGO_TARGET_DIR="${target_dir}" cargo run --quiet --bin gatk-rs -- \
  HaplotypeCaller \
  -R parity/fixtures/reference.fa \
  -I parity/fixtures/sample.bam \
  -O "${tmp_vcf}" \
  -L chr1:1-32

python3 - "${tmp_vcf}" "${json_out}" <<'PY'
import json
import pathlib
import sys

vcf = pathlib.Path(sys.argv[1])
out = pathlib.Path(sys.argv[2])

if not vcf.exists():
    print("[p11-activation] missing VCF output", file=sys.stderr)
    raise SystemExit(1)

lines = vcf.read_text(encoding="utf-8", errors="replace").splitlines()
has_header = any(l.startswith("##fileformat=") for l in lines) and any(
    l.startswith("#CHROM") for l in lines
)
variants = [l for l in lines if l and not l.startswith("#")]

# Pass when the VCF is structurally valid; variant rows optional on small fixtures.
status = "pass" if has_header else "fail"
payload = {
    "label": "phase11-hc-output-activation-contracts",
    "status": status,
    "has_vcf_header": bool(has_header),
    "variant_record_count": len(variants),
    "notes": (
        "Default assembly-region-v1 pipeline emitted variant record(s)"
        if len(variants) > 0
        else "Valid assembly-region-v1 header; zero variant rows in fixture window"
    ),
}
out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
if not has_header:
    print("[p11-activation] invalid VCF header", file=sys.stderr)
    raise SystemExit(1)
print(f"[p11-activation] status={status} variant_record_count={len(variants)}")
PY

rm -f "${tmp_vcf}"
