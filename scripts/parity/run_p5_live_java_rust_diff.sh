#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
report_dir="${repo_root}/parity/reports"
manifest="${P5_LIVE_MANIFEST:-${repo_root}/parity/fixtures/p5_live_regions.tsv}"
mkdir -p "${report_dir}"

gatk_image="${GATK_DOCKER_IMAGE:-us.gcr.io/broad-gatk/gatk:4.4.0.0}"
gatk_platform="${GATK_DOCKER_PLATFORM:-linux/amd64}"

passed=0
failed=0
rows_json="${report_dir}/p5_live_java_rust_diff_rows.jsonl"
: > "${rows_json}"

while IFS=$'\t' read -r case_id reference_fa sam_fixture interval; do
  [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
  ref_path="${repo_root}/parity/fixtures/${reference_fa}"
  sam_path="${repo_root}/parity/fixtures/${sam_fixture}"
  cp "${ref_path}" "${report_dir}/${case_id}.ref.fa"
  samtools faidx "${report_dir}/${case_id}.ref.fa"
  rm -f "${report_dir}/${case_id}.ref.dict"
  docker run --rm --platform "${gatk_platform}" -v "${repo_root}:/app" -w /app "${gatk_image}" \
    bash -c "gatk CreateSequenceDictionary -R '/app/parity/reports/${case_id}.ref.fa' -O '/app/parity/reports/${case_id}.ref.dict'" \
    >/dev/null 2>&1

  tmp_bam="${report_dir}/${case_id}.live.bam"
  java_out="${report_dir}/${case_id}.java.live.out"
  vcf_out="${report_dir}/${case_id}.java.live.vcf"
  tmp_bam_in_container="/app/parity/reports/${case_id}.live.bam"
  vcf_out_in_container="/app/parity/reports/${case_id}.java.live.vcf"
  samtools view -bS "${sam_path}" > "${tmp_bam}"
  samtools index "${tmp_bam}"

  set +e
  docker run --rm --platform "${gatk_platform}" -v "${repo_root}:/app" -w /app "${gatk_image}" \
    bash -c "gatk HaplotypeCaller -R '/app/parity/reports/${case_id}.ref.fa' -I '${tmp_bam_in_container}' -O '${vcf_out_in_container}' -L '${interval}' --debug-assembly true" \
    > "${java_out}" 2>&1
  java_code=$?
  set -e
  if [[ "${java_code}" -ne 0 ]]; then
    failed=$((failed + 1))
    python3 - <<PY >> "${rows_json}"
import json
print(json.dumps({"case":"${case_id}","interval":"${interval}","ok":False,"reason":"java_exit_${java_code}"}))
PY
    continue
  fi

  set +e
  P5_LIVE_SAM="${sam_path}" P5_LIVE_JAVA_OUT="${java_out}" P5_LIVE_JAVA_VCF="${vcf_out}" \
    cargo test -p gatk-haplotypecaller --test p5_live_java_rust_diff_test --locked live_java_eventmap_haplotype_signatures_cover_rust_candidates >/dev/null 2>&1
  rust_code=$?
  set -e
  if [[ "${rust_code}" -eq 0 ]]; then
    passed=$((passed + 1))
    python3 - <<PY >> "${rows_json}"
import json
print(json.dumps({"case":"${case_id}","interval":"${interval}","ok":True}))
PY
  else
    failed=$((failed + 1))
    python3 - <<PY >> "${rows_json}"
import json
print(json.dumps({"case":"${case_id}","interval":"${interval}","ok":False,"reason":"rust_compare_failed"}))
PY
  fi
done < "${manifest}"

summary="${report_dir}/p5_live_java_rust_diff_summary.json"
python3 - <<PY
import json, pathlib
rows=[json.loads(l) for l in pathlib.Path("${rows_json}").read_text(encoding="utf-8").splitlines() if l.strip()]
passed=sum(1 for r in rows if r.get("ok"))
failed=sum(1 for r in rows if not r.get("ok"))
summary={
  "label":"phase5-live-java-rust-diff",
  "total":len(rows),
  "passed":passed,
  "failed":failed,
  "pass_rate": (passed/len(rows)) if rows else 0.0,
  "rows":rows
}
pathlib.Path("${summary}").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(f"[p5-live-diff] wrote ${summary}")
print(f"[p5-live-diff] passed={passed} failed={failed}")
if failed:
  raise SystemExit(1)
PY
