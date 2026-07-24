#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
report_dir="${repo_root}/parity/reports"
mkdir -p "${report_dir}"
cd "${repo_root}"

log_file="${report_dir}/p14_multidataset_equivalence.log"

log() {
  printf '[%s] %s\n' "$(date -u +"%Y-%m-%dT%H:%M:%SZ")" "$1" | tee -a "${log_file}"
}

run_cmd() {
  log "RUN: $*"
  "$@" 2>&1 | tee -a "${log_file}"
}

write_pending_case() {
  local case_id="$1"
  local reason="$2"
  local env_hint="$3"
  python3 - "${report_dir}/p14_${case_id}.json" "${report_dir}/p14_${case_id}.md" "${case_id}" "${reason}" "${env_hint}" <<'PY'
import json
import pathlib
import sys

json_out = pathlib.Path(sys.argv[1])
md_out = pathlib.Path(sys.argv[2])
case_id = sys.argv[3]
reason = sys.argv[4]
env_hint = sys.argv[5]

payload = {
    "label": "phase14-multidataset-equivalence-case",
    "case_id": case_id,
    "status": "pending_data",
    "reason": reason,
    "env_hint": env_hint,
}
json_out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
md_out.write_text(
    "\n".join(
        [
            f"# P14 Case {case_id}",
            "",
            "- status: **pending_data**",
            f"- reason: `{reason}`",
            f"- action: `{env_hint}`",
        ]
    )
    + "\n",
    encoding="utf-8",
)
PY
}

run_generic_case() {
  local case_id="$1"
  local reference="$2"
  local bam="$3"
  local interval="$4"
  local truth_vcf="$5"
  local truth_bed="$6"
  local chrom="$7"
  local java_vcf="${report_dir}/p14_${case_id}.java.vcf"
  local rust_vcf="${report_dir}/p14_${case_id}.rust.vcf"
  local target_dir="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target-parity}"
  mkdir -p "${target_dir}"

  set +e
  docker run --rm --platform "${GATK_DOCKER_PLATFORM:-linux/amd64}" \
    -v "${repo_root}:${repo_root}" \
    -w "${repo_root}" \
    "${GATK_DOCKER_IMAGE:-us.gcr.io/broad-gatk/gatk:4.4.0.0}" \
    gatk HaplotypeCaller \
    -R "${reference}" \
    -I "${bam}" \
    -O "${java_vcf}" \
    -L "${interval}" \
    --verbosity ERROR >/dev/null 2>&1
  local java_exit=$?

  CARGO_TARGET_DIR="${target_dir}" cargo run --quiet --bin gatk-rs -- \
    HaplotypeCaller \
    -R "${reference}" \
    -I "${bam}" \
    -O "${rust_vcf}" \
    -L "${interval}" >/dev/null 2>&1
  local rust_exit=$?
  set -e

  P13_JAVA_VCF="${java_vcf}" \
  P13_RUST_VCF="${rust_vcf}" \
  P13_TRUTH_VCF="${truth_vcf}" \
  P13_REGIONS_BED="${truth_bed}" \
  P13_CHROM="${chrom}" \
  P13_EVAL_INTERVAL="${interval}" \
  run_cmd ./scripts/parity/run_p13_realworld_truth_eval.sh

  python3 - \
    "${report_dir}/p14_${case_id}.json" \
    "${report_dir}/p14_${case_id}.md" \
    "${report_dir}/p13_realworld_truth_eval.json" \
    "${case_id}" \
    "${java_exit}" \
    "${rust_exit}" \
    <<'PY'
import json
import pathlib
import sys

json_out = pathlib.Path(sys.argv[1])
md_out = pathlib.Path(sys.argv[2])
p13_json = pathlib.Path(sys.argv[3])
case_id = sys.argv[4]
java_exit = int(sys.argv[5])
rust_exit = int(sys.argv[6])

p13 = json.loads(p13_json.read_text(encoding="utf-8"))
payload = {
    "label": "phase14-multidataset-equivalence-case",
    "case_id": case_id,
    "status": "pass" if java_exit == 0 and rust_exit == 0 and p13.get("status") == "pass" else "needs_attention",
    "java_exit": java_exit,
    "rust_exit": rust_exit,
    "truth_eval": p13,
}
json_out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
md_out.write_text(
    "\n".join(
        [
            f"# P14 Case {case_id}",
            "",
            f"- status: **{payload['status']}**",
            f"- java/rust exit: `{java_exit}/{rust_exit}`",
            f"- truth eval status: `{p13.get('status')}`",
            f"- java f1: `{(p13.get('java') or {}).get('f1', 0):.4f}`",
            f"- rust f1: `{(p13.get('rust') or {}).get('f1', 0):.4f}`",
        ]
    )
    + "\n",
    encoding="utf-8",
)
PY
}

log "=== Phase 14 multi-dataset equivalence start ==="

# Case 1: NA12878 + GIAB (always runnable via existing P12/P13 full harness)
run_cmd ./scripts/parity/run_p12_p13_realworld_full.sh
python3 - "${report_dir}/p14_na12878_20k_giab_b37.json" "${report_dir}/p14_na12878_20k_giab_b37.md" "${report_dir}/p12_p13_realworld_full_summary.json" <<'PY'
import json
import pathlib
import sys

json_out = pathlib.Path(sys.argv[1])
md_out = pathlib.Path(sys.argv[2])
src = pathlib.Path(sys.argv[3])
payload = json.loads(src.read_text(encoding="utf-8"))
case_payload = {
    "label": "phase14-multidataset-equivalence-case",
    "case_id": "na12878_20k_giab_b37",
    "status": payload.get("status", "needs_attention"),
    "source_summary": payload,
}
json_out.write_text(json.dumps(case_payload, indent=2) + "\n", encoding="utf-8")
md_out.write_text(
    "\n".join(
        [
            "# P14 Case na12878_20k_giab_b37",
            "",
            f"- status: **{case_payload['status']}**",
            f"- p12 status: `{payload.get('p12_status')}`",
            f"- p13 status: `{payload.get('p13_status')}`",
            f"- java/rust variants: `{(payload.get('p12') or {}).get('java_variant_count')} / {(payload.get('p12') or {}).get('rust_variant_count')}`",
        ]
    )
    + "\n",
    encoding="utf-8",
)
PY

# Case 2: Syndip CHM synthetic diploid (optional, requires user-provided BAM/REF).
if [[ -n "${P14_SYNDIP_REFERENCE:-}" && -n "${P14_SYNDIP_BAM:-}" && -n "${P14_SYNDIP_TRUTH_VCF:-}" && -n "${P14_SYNDIP_TRUTH_BED:-}" && -n "${P14_SYNDIP_INTERVAL:-}" ]]; then
  run_generic_case \
    "syndip_chm37" \
    "${P14_SYNDIP_REFERENCE}" \
    "${P14_SYNDIP_BAM}" \
    "${P14_SYNDIP_INTERVAL}" \
    "${P14_SYNDIP_TRUTH_VCF}" \
    "${P14_SYNDIP_TRUTH_BED}" \
    "${P14_SYNDIP_CHROM:-20}"
else
  write_pending_case \
    "syndip_chm37" \
    "missing required Syndip environment variables" \
    "set P14_SYNDIP_REFERENCE,P14_SYNDIP_BAM,P14_SYNDIP_TRUTH_VCF,P14_SYNDIP_TRUTH_BED,P14_SYNDIP_INTERVAL"
fi

# Case 3: precisionFDA HG003 (optional, requires user-provided assets).
if [[ -n "${P14_PFDA_HG003_REFERENCE:-}" && -n "${P14_PFDA_HG003_BAM:-}" && -n "${P14_PFDA_HG003_TRUTH_VCF:-}" && -n "${P14_PFDA_HG003_TRUTH_BED:-}" && -n "${P14_PFDA_HG003_INTERVAL:-}" ]]; then
  run_generic_case \
    "precisionfda_hg003_grch38" \
    "${P14_PFDA_HG003_REFERENCE}" \
    "${P14_PFDA_HG003_BAM}" \
    "${P14_PFDA_HG003_INTERVAL}" \
    "${P14_PFDA_HG003_TRUTH_VCF}" \
    "${P14_PFDA_HG003_TRUTH_BED}" \
    "${P14_PFDA_HG003_CHROM:-20}"
else
  write_pending_case \
    "precisionfda_hg003_grch38" \
    "missing required precisionFDA HG003 environment variables" \
    "set P14_PFDA_HG003_REFERENCE,P14_PFDA_HG003_BAM,P14_PFDA_HG003_TRUTH_VCF,P14_PFDA_HG003_TRUTH_BED,P14_PFDA_HG003_INTERVAL"
fi

# Case 4: precisionFDA HG004 (optional, requires user-provided assets).
if [[ -n "${P14_PFDA_HG004_REFERENCE:-}" && -n "${P14_PFDA_HG004_BAM:-}" && -n "${P14_PFDA_HG004_TRUTH_VCF:-}" && -n "${P14_PFDA_HG004_TRUTH_BED:-}" && -n "${P14_PFDA_HG004_INTERVAL:-}" ]]; then
  run_generic_case \
    "precisionfda_hg004_grch38" \
    "${P14_PFDA_HG004_REFERENCE}" \
    "${P14_PFDA_HG004_BAM}" \
    "${P14_PFDA_HG004_INTERVAL}" \
    "${P14_PFDA_HG004_TRUTH_VCF}" \
    "${P14_PFDA_HG004_TRUTH_BED}" \
    "${P14_PFDA_HG004_CHROM:-20}"
else
  write_pending_case \
    "precisionfda_hg004_grch38" \
    "missing required precisionFDA HG004 environment variables" \
    "set P14_PFDA_HG004_REFERENCE,P14_PFDA_HG004_BAM,P14_PFDA_HG004_TRUTH_VCF,P14_PFDA_HG004_TRUTH_BED,P14_PFDA_HG004_INTERVAL"
fi

run_cmd python3 ./scripts/parity/build_p14_equivalence_summary.py
log "=== Phase 14 multi-dataset equivalence complete ==="
