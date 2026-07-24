#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
report_dir="${repo_root}/parity/reports"
mkdir -p "${report_dir}"
cd "${repo_root}"

java_vcf="${P13_JAVA_VCF:-${report_dir}/p12_realworld_na12878_20k.java.vcf}"
rust_vcf="${P13_RUST_VCF:-${report_dir}/p12_realworld_na12878_20k.rust.vcf}"
truth_vcf="${P13_TRUTH_VCF:-}"
chrom_filter="${P13_CHROM:-}"
regions_bed="${P13_REGIONS_BED:-}"
thresholds_json="${P13_THRESHOLDS_JSON:-}"
json_out="${report_dir}/p13_realworld_truth_eval.json"
md_out="${report_dir}/p13_realworld_truth_eval.md"
strict_gate="${P13_STRICT_GATE:-0}"

if [[ -z "${truth_vcf}" ]]; then
  python3 - "${json_out}" "${md_out}" <<'PY'
import json
import pathlib
import sys

json_out = pathlib.Path(sys.argv[1])
md_out = pathlib.Path(sys.argv[2])
payload = {
    "label": "phase13-realworld-truth-eval",
    "status": "truth_missing",
    "notes": "Set P13_TRUTH_VCF to a GIAB truth VCF to evaluate Java and Rust callsets against external truth.",
}
json_out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
md_out.write_text(
    "# P13 Real-world Truth Eval\n\n- status: **truth_missing**\n- next: export `P13_TRUTH_VCF=/path/to/giab_truth.vcf.gz` and rerun\n",
    encoding="utf-8",
)
print("[p13-truth] truth VCF missing; evaluation skipped")
PY
  exit 0
fi

# Restrict truth + calls to the HC interval when set (slice runs). Prefer P13_EVAL_INTERVAL, else P12_INTERVAL.
eval_interval="${P13_EVAL_INTERVAL:-${P12_INTERVAL:-}}"
regions_arg=()
if [[ -n "${regions_bed}" ]]; then
  regions_arg=(--regions-bed "${regions_bed}")
fi

python3 "${repo_root}/scripts/parity/p13_truth_eval.py" \
  --java-vcf "${java_vcf}" \
  --rust-vcf "${rust_vcf}" \
  --truth-vcf "${truth_vcf}" \
  --json-out "${json_out}" \
  --md-out "${md_out}" \
  --chrom-filter "${chrom_filter}" \
  $(if [[ ${#regions_arg[@]} -gt 0 ]]; then printf '%s ' "${regions_arg[@]}"; fi) \
  --eval-interval "${eval_interval}" \
  $(if [[ -n "${thresholds_json}" && -f "${thresholds_json}" ]]; then echo --thresholds-json "${thresholds_json}"; fi) \
  $(if [[ "${strict_gate}" == "1" ]]; then echo --strict-gate; fi)
