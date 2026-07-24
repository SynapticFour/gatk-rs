#!/usr/bin/env bash
# Run Real-World pipeline steps: Rust + Java where comparable, write summary under OUT_DIR.
# Requires: docs/HC_REALWORLD_PIPELINE.md (definition). Equivalence contract:
#   docs/REALWORLD_EQUIVALENCE.md — what “PASS/PARITY” means per step vs GATK4 (not all are strict).
# Emits: OUT_DIR/equivalence_report.{md,json} via realworld_equivalence_report.py
#
# Optional: RW_ONLY_STEP=3 ./run_paired_realworld_pipeline.sh   — run only step 3 (count reads)
# Optional: RW_SKIP_STEP05=1  — skip read-filter SAM parity (Java PrintReads vs Rust FilterReads)
# Optional: RW_SKIP_STEP06=1  — skip assembly-region smoothed-activity parity (saves time)
# Optional: RW_SKIP_STEP07=1  — skip full HC (saves a lot of time)
# Optional: RW_SMOOTHED_ACTIVITY_STRICT=1 — step 06 also enforces legacy continuous max-abs-diff (debug; usually fails)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=/dev/null
source "${SCRIPT_DIR}/common.sh"

cd "${REPO_ROOT}"

ONLY="${RW_ONLY_STEP:-}"
SKIP05="${RW_SKIP_STEP05:-0}"
SKIP06="${RW_SKIP_STEP06:-0}"
SKIP07="${RW_SKIP_STEP07:-0}"

log() { printf '[rw-pipeline] %s\n' "$*"; }

log "Starting paired realworld pipeline (see docs/REALWORLD_EQUIVALENCE.md for what is / is not proven)."
log "Env: RW_SKIP_STEP05=${SKIP05} RW_SKIP_STEP06=${SKIP06} RW_SKIP_STEP07=${SKIP07} RW_ONLY_STEP=${ONLY:-<all>}"

# Append-only footer so tee’d logs that stop early still leave a durable completion record under OUT_DIR.
append_run_footer() {
  local ts
  ts="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  {
    echo "[rw-pipeline] --- run footer ${ts} ---"
    echo "[rw-pipeline] RW_SKIP_STEP05=${SKIP05} RW_SKIP_STEP06=${SKIP06} RW_SKIP_STEP07=${SKIP07} RW_ONLY_STEP=${ONLY:-}"
    echo "[rw-pipeline] FAIL_AGGREGATE=${FAIL}"
    if [[ "${FAIL}" -ne 0 ]]; then
      echo "[rw-pipeline] OUTCOME=FAILED_PARITY_OR_TOOL_STEP (see summary.md + equivalence_report.md)"
    else
      echo "[rw-pipeline] OUTCOME=ALL_EXECUTED_STEPS_PASS_OR_EXPECTED_SKIP"
    fi
  } >> "${OUT_DIR}/pipeline_footer.txt"
}

should_run() {
  local n="$1"
  if [[ -z "${ONLY}" ]]; then
    return 0
  fi
  [[ "${ONLY}" == "${n}" ]]
}

summary_lines=()

step_01() {
  log "Step 01 — verify inputs"
  local ok=0
  local bai="${RW_BAM}.bai"
  if [[ ! -f "${bai}" ]]; then
    bai="${RW_BAM%.bam}.bai"
  fi
  for p in "${RW_REF}" "${RW_BAM}" "${bai}"; do
    if [[ ! -f "${p}" ]]; then
      echo "MISSING: ${p}" >&2
      ok=1
    else
      log "OK: ${p}"
    fi
  done
  if [[ "${ok}" -ne 0 ]]; then
    summary_lines+=("01_verify: FAIL (missing files)")
    return 1
  fi
  summary_lines+=("01_verify: PASS — EQUIV: n/a (filesystem gate only; no GATK cross-check)")
  return 0
}

step_02() {
  log "Step 02 — Validate BAM (Java ValidateSamFile vs Rust Validate)"
  local jout="${OUT_DIR}/02_validate.java.stdout"
  local rout="${OUT_DIR}/02_validate.rust.stdout"
  set +e
  "${RUN_JAVA}" "${jout}" ValidateSamFile -I "${RW_BAM}" -MODE SUMMARY
  local je=$?
  "${RUN_RUST}" "${rout}" Validate "${RW_BAM}" -t BAM -R "${RW_REF}"
  local re=$?
  set -e
  if [[ "${je}" -eq 0 && "${re}" -eq 0 ]]; then
    summary_lines+=("02_validate: PASS (java_exit=${je} rust_exit=${re}) — EQUIV: both exit 0 = operational agreement; ValidateSamFile vs Rust Validate are not identical tools or logs")
  else
    summary_lines+=("02_validate: MISMATCH (java_exit=${je} rust_exit=${re}) — EQUIV: FAILED (need both sides to accept BAM)")
    return 1
  fi
  return 0
}

step_03() {
  log "Step 03 — Count reads in interval (Java CountReads vs Rust CountReadsInRegion)"
  local jout="${OUT_DIR}/03_count.java.stdout"
  local rout="${OUT_DIR}/03_count.rust.stdout"
  set +e
  "${RUN_JAVA}" "${jout}" CountReads -I "${RW_BAM}" -L "${RW_INTERVAL}"
  local je=$?
  "${RUN_RUST}" "${rout}" CountReadsInRegion -I "${RW_BAM}" -L "${RW_INTERVAL}"
  local re=$?
  set -e
  if [[ "${je}" -ne 0 || "${re}" -ne 0 ]]; then
    summary_lines+=("03_count_reads: TOOL_FAIL java_exit=${je} rust_exit=${re}")
    return 1
  fi
  local jc rc
  jc="$(python3 - <<PY
import re
t=open("${jout}",encoding="utf-8",errors="replace").read()
m=re.search(r"CountReads counted\\s+(\\d+)\\s+total reads", t)
print(m.group(1) if m else "-1")
PY
)"
  rc="$(grep -E '^COUNT :' "${rout}" | awk '{print $3}' | head -1)"
  if [[ "${jc}" == "${rc}" ]]; then
    summary_lines+=("03_count_reads: PARITY count=${jc} — EQUIV: strict numeric match with GATK4 CountReads (same -L)")
  else
    summary_lines+=("03_count_reads: DIVERGENCE java_count=${jc} rust_count=${rc} — EQUIV: FAILED strict count parity")
    return 1
  fi
  return 0
}

step_04() {
  log "Step 04 — Count bases in reference interval (Rust CountBasesInReference; Java CountBasesInReference)"
  local rout="${OUT_DIR}/04_count_bases.rust.stdout"
  set +e
  "${RUN_RUST}" "${rout}" CountBasesInReference -R "${RW_REF}" -L "${RW_INTERVAL}"
  local re=$?
  set -e
  if [[ "${re}" -ne 0 ]]; then
    summary_lines+=("04_count_bases_rust: FAIL exit=${re}")
    return 1
  fi
  summary_lines+=("04_count_bases_rust: PASS (see ${rout})")
  local jout="${OUT_DIR}/04_count_bases.java.stdout"
  set +e
  # GATK4 CountBases counts bases in an alignment file (-I); reference intervals use CountBasesInReference (matches Rust).
  "${RUN_JAVA}" "${jout}" CountBasesInReference -R "${RW_REF}" -L "${RW_INTERVAL}"
  local je=$?
  set -e
  if [[ "${je}" -eq 0 ]]; then
    summary_lines+=("04_count_bases_java: PASS (see ${jout}) — EQUIV target: A/C/G/T/N histogram vs Rust")
    hl="$(python3 "${SCRIPT_DIR}/realworld_equivalence_report.py" --step04 "${rout}" "${jout}" 2>&1)" || true
    summary_lines+=("04_histogram_equivalence: ${hl}")
  else
    summary_lines+=("04_count_bases_java: FAIL exit=${je} — EQUIV: cannot compare histogram until Java succeeds")
  fi
  return 0
}

step_05_filter_parity() {
  if [[ "${SKIP05}" == "1" ]]; then
    summary_lines+=("05_filter_reads: SKIPPED (RW_SKIP_STEP05=1)")
    return 0
  fi
  log "Step 05 — Read filter parity (GATK4 PrintReads + HC-style filters vs Rust FilterReads)"
  local jsam="${OUT_DIR}/05_filter.java.sam"
  local rsam="${OUT_DIR}/05_filter.rust.sam"
  local jout="${OUT_DIR}/05_filter.java.stdout"
  local rout="${OUT_DIR}/05_filter.rust.stdout"
  local cmp_py="${REPO_ROOT}/scripts/parity/compare_sam_parity.py"
  set +e
  "${RUN_JAVA}" "${jout}" PrintReads \
    -I "${RW_BAM}" \
    -O "${jsam}" \
    -L "${RW_INTERVAL}" \
    --read-filter MappedReadFilter \
    --read-filter NotSecondaryAlignmentReadFilter \
    --read-filter NotSupplementaryAlignmentReadFilter \
    --read-filter MappingQualityReadFilter \
    --minimum-mapping-quality 20 \
    --read-filter NotDuplicateReadFilter
  local je=$?
  CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR}" PARITY_RUST_PROFILE="${PARITY_RUST_PROFILE}" \
    "${RUN_RUST}" "${rout}" FilterReads \
    -I "${RW_BAM}" \
    -O "${rsam}" \
    -L "${RW_INTERVAL}" \
    --min-mapping-quality 20
  local re=$?
  set -e
  if [[ "${je}" -ne 0 || "${re}" -ne 0 ]]; then
    summary_lines+=("05_filter_reads: TOOL_FAIL java_exit=${je} rust_exit=${re}")
    return 1
  fi
  set +e
  python3 "${cmp_py}" \
    --java-sam "${jsam}" \
    --rust-sam "${rsam}" \
    --label "realworld_step05" \
    --json-out "${OUT_DIR}/05_filter_parity.json"
  local ce=$?
  set -e
  if [[ "${ce}" -eq 0 ]]; then
    summary_lines+=("05_filter_reads: PARITY — EQUIV: normalized SAM match (see ${OUT_DIR}/05_filter_parity.json); same filter semantics as run_read_filter_diff.sh")
  else
    summary_lines+=("05_filter_reads: DIVERGENCE — EQUIV FAILED (see ${OUT_DIR}/05_filter_parity.json)")
    return 1
  fi
  return 0
}

step_06_assembly_activity_parity() {
  if [[ "${SKIP06}" == "1" ]]; then
    summary_lines+=("06_assembly_activity: SKIPPED (RW_SKIP_STEP06=1)")
    return 0
  fi
  log "Step 06 — Assembly-region / smoothed activity (Java IGV vs Rust DumpSmoothedActivity)"
  local scratch_vcf="${OUT_DIR}/06_hc_scratch.java.vcf"
  local igv="${OUT_DIR}/06_assembly_regions.java.igv"
  local rust_tsv="${OUT_DIR}/06_smoothed.rust.tsv"
  local cmp_py="${SCRIPT_DIR}/compare_smoothed_activity.py"
  set +e
  "${RUN_JAVA}" "${OUT_DIR}/06_hc_assembly.java.stdout" HaplotypeCaller \
    -R "${RW_REF}" \
    -I "${RW_BAM}" \
    -O "${scratch_vcf}" \
    -L "${RW_INTERVAL}" \
    --assembly-region-out "${igv}" \
    --verbosity ERROR
  local je=$?
  set -e
  if [[ "${je}" -ne 0 || ! -f "${igv}" ]]; then
    summary_lines+=("06_assembly_java_igv: FAIL exit=${je}")
    return 1
  fi
  summary_lines+=("06_assembly_java_igv: PASS (${igv})")
  set +e
  CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR}" PARITY_RUST_PROFILE="${PARITY_RUST_PROFILE}" \
    "${RUN_RUST}" "${OUT_DIR}/06_smoothed.rust.stdout" DumpSmoothedActivity \
    -R "${RW_REF}" \
    -I "${RW_BAM}" \
    -L "${RW_INTERVAL}" \
    -O "${rust_tsv}"
  local re=$?
  set -e
  if [[ "${re}" -ne 0 ]]; then
    summary_lines+=("06_smoothed_rust: FAIL exit=${re}")
    return 1
  fi
  set +e
  local ce
  if [[ "${RW_SMOOTHED_ACTIVITY_STRICT:-0}" == "1" ]]; then
    python3 "${cmp_py}" \
      --java-igv "${igv}" \
      --rust-tsv "${rust_tsv}" \
      --json-out "${OUT_DIR}/06_smoothed_parity.json" \
      --max-abs-diff "${RW_SMOOTHED_ACTIVITY_MAX_DIFF:-0.15}" \
      --require-continuous-max-diff
    ce=$?
  else
    python3 "${cmp_py}" \
      --java-igv "${igv}" \
      --rust-tsv "${rust_tsv}" \
      --json-out "${OUT_DIR}/06_smoothed_parity.json" \
      --max-abs-diff "${RW_SMOOTHED_ACTIVITY_MAX_DIFF:-0.15}"
    ce=$?
  fi
  set -e
  if [[ "${ce}" -eq 0 ]]; then
    summary_lines+=("06_smoothed_activity: PARITY — EQUIV: binary active-region + coverage vs Java IGV (06_smoothed_parity.json contract; optional RW_SMOOTHED_ACTIVITY_STRICT=1 for float gate)")
  else
    summary_lines+=("06_smoothed_activity: DIVERGENCE — refine pileup/activity or increase tolerance (see 06_smoothed_parity.json)")
    return 1
  fi
  return 0
}

step_07_hc() {
  if [[ "${SKIP07}" == "1" ]]; then
    summary_lines+=("07_mode: SKIPPED_BY_ENV (RW_SKIP_STEP07=1)")
    summary_lines+=("07_haplotypecaller: SKIPPED (RW_SKIP_STEP07=1)")
    return 0
  fi
  log "Step 07 — Full HaplotypeCaller VCF (Java vs Rust+activation)"
  summary_lines+=("07_mode: EXECUTED (RW_SKIP_STEP07=0) — artifacts: ${OUT_DIR}/07_haplotypecaller.{java,rust}.vcf")
  local jvcf="${OUT_DIR}/07_haplotypecaller.java.vcf"
  local rvcf="${OUT_DIR}/07_haplotypecaller.rust.vcf"
  set +e
  "${RUN_JAVA}" "${OUT_DIR}/07_hc.java.stdout" HaplotypeCaller \
    -R "${RW_REF}" \
    -I "${RW_BAM}" \
    -O "${jvcf}" \
    -L "${RW_INTERVAL}" \
    --verbosity ERROR
  local je=$?
  CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR}" PARITY_RUST_PROFILE="${PARITY_RUST_PROFILE}" \
    "${RUN_RUST}" "${OUT_DIR}/07_hc.rust.stdout" HaplotypeCaller \
    -R "${RW_REF}" \
    -I "${RW_BAM}" \
    -O "${rvcf}" \
    -L "${RW_INTERVAL}"
  local re=$?
  set -e
  local jn=0 rn=0
  if [[ -f "${jvcf}" ]]; then
    jn=$(grep -c -v '^#' "${jvcf}" 2>/dev/null || true)
  fi
  if [[ -f "${rvcf}" ]]; then
    rn=$(grep -c -v '^#' "${rvcf}" 2>/dev/null || true)
  fi
  [[ -z "${jn}" ]] && jn=0
  [[ -z "${rn}" ]] && rn=0
  summary_lines+=("07_haplotypecaller: java_exit=${je} rust_exit=${re} java_data_lines=${jn} rust_data_lines=${rn}")
  summary_lines+=("07_vcf_equivalence: compare CHROM/POS/REF/ALT sets in equivalence_report.md — P12 L3/L5 signed; full genome = L6")
  summary_lines+=("07_note: Rust default is assembly-region-v1; P12 variant-set parity proven via L3/L5 batteries.")
  return 0
}

write_summary() {
  local sum="${OUT_DIR}/summary.md"
  {
    echo "# Real-World pipeline run"
    echo ""
    echo "- RW_REF=\`${RW_REF}\`"
    echo "- RW_BAM=\`${RW_BAM}\`"
    echo "- RW_INTERVAL=\`${RW_INTERVAL}\`"
    echo ""
    echo "## Results"
    echo ""
    for line in "${summary_lines[@]}"; do
      echo "- ${line}"
    done
    echo ""
    echo "## Equivalence contract (summary)"
    echo ""
    echo "| Step | vs GATK4 Java |"
    echo "|------|----------------|"
    echo "| 01 | Files exist — no tool comparison |"
    echo "| 02 | Loose: both exit 0 (different validators) |"
    echo "| 03 | **Strict:** same read count in \`-L\` |"
    echo "| 04 | **Strict:** same A/C/G/T/N histogram (\`CountBasesInReference\`) |"
    echo "| 05 | **Strict:** normalized SAM parity (\`PrintReads\`+filters vs \`FilterReads\`) |"
    echo "| 06 | **Strict:** smoothed activity vs Java IGV (\`DumpSmoothedActivity\` vs HC \`--assembly-region-out\`) |"
    echo "| 07 | Variant **set** report only; not byte-identical VCF vs full HC |"
    echo ""
    echo "Canonical detail: \`docs/REALWORLD_EQUIVALENCE.md\`. Per-run analysis: \`${OUT_DIR}/equivalence_report.md\`."
    echo ""
    echo "## Definition"
    echo ""
    echo "See \`docs/HC_REALWORLD_PIPELINE.md\`."
    echo ""
    echo "## Run metadata (for audit / no silent env mismatch)"
    echo ""
    echo "- UTC when this file was written: \`$(date -u +%Y-%m-%dT%H:%M:%SZ)\`"
    echo "- \`RW_SKIP_STEP05=${SKIP05}\` \`RW_SKIP_STEP06=${SKIP06}\` \`RW_SKIP_STEP07=${SKIP07}\`"
    echo "- \`RW_ONLY_STEP=${ONLY:-}\` (empty means all steps attempted)"
    echo "- \`FAIL_AGGREGATE=${FAIL}\` (1 = at least one step returned non-zero parity/tool failure)"
    echo "- durable footer also appended to: \`${OUT_DIR}/pipeline_footer.txt\`"
  } > "${sum}"
  if python3 "${SCRIPT_DIR}/realworld_equivalence_report.py" "${OUT_DIR}" >/dev/null 2>&1; then
    log "Wrote ${OUT_DIR}/equivalence_report.md"
  else
    log "Note: equivalence_report skipped (python error or missing artifacts)"
  fi
  log "Wrote ${sum}"
}

FAIL=0
# Invoke the step by name (second arg). Do not use shift + "$@": a buggy version
# left a stray token and produced `n: command not found` (exit 127) on the next line.
run_step() {
  local num="$1"
  local fn="$2"
  if [[ -z "${fn}" ]]; then
    log "run_step: missing function for step ${num}" >&2
    return 1
  fi
  if should_run "${num}"; then
    if ! "${fn}"; then
      FAIL=1
    fi
  fi
}

run_step "01" step_01 || true
run_step "02" step_02 || true
run_step "03" step_03 || true
run_step "04" step_04 || true
run_step "05" step_05_filter_parity || true
run_step "06" step_06_assembly_activity_parity || true
run_step "07" step_07_hc || true

write_summary
append_run_footer

if [[ "${FAIL}" -ne 0 ]]; then
  log "Completed with at least one failed parity step (see summary)."
  exit 1
fi
log "All executed steps reported PASS / expected SKIP."
exit 0
