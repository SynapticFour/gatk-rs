#!/usr/bin/env bash
# L2 gate: fresh Rust dump vs frozen Java dump (parity/fixtures/hc-full-parity/java_dumps/).
# Phase C paths: c1 raw-activity, c2 smoothed, c3 active-locus, c4-gl, c5-multi (multisample raw-activity).
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=lib_pinned_gatk.sh
source "${script_dir}/lib_pinned_gatk.sh"

repo_root="${GATK_RS_REPO_ROOT}"
cd "${repo_root}"
export CARGO_TARGET_DIR="${PARITY_CARGO_TARGET_DIR:-${repo_root}/target}"

# Production L2 gate: Phase-E harness flags must not leak from the developer shell.
unset P12_PHASE_E || true
unset P12_BASELINE_EMIT_FILTER || true

pin_short="${GATK_PINNED_SHA_SHORT:-2dbc0258}"
java_root="${repo_root}/parity/fixtures/hc-full-parity/java_dumps"
report_dir="${repo_root}/parity/reports/hc-full-parity-l2"
tmp_dir="${report_dir}/tmp"
mkdir -p "${tmp_dir}" "${report_dir}"

compare_py="${script_dir}/compare_hc_full_parity_l2.py"
dump_java="${script_dir}/run_hc_full_parity_java_dump.sh"
strict="${PARITY_HC_FULL_L2_STRICT:-1}"

resolve_alignment_path_l2() {
  local path="$1"
  if [[ "${path}" != *.sam ]]; then
    printf '%s\n' "${path}"
    return
  fi
  local cache_dir="${repo_root}/parity/build/sam-indexed-bam"
  mkdir -p "${cache_dir}"
  local out="${cache_dir}/$(basename "${path%.sam}").bam"
  if [[ ! -f "${out}" ]]; then
    echo "[hc-full-parity-l2] sam->bam $(basename "${path}") ..." >&2
    local sam_err="${cache_dir}/$(basename "${path%.sam}").samtools.err"
    if ! samtools view -bS "${path}" 2>"${sam_err}" | samtools sort -o "${out}" 2>>"${sam_err}"; then
      echo "[hc-full-parity-l2] sam->bam failed for $(basename "${path}") (see ${sam_err})" >&2
      return 1
    fi
    echo "[hc-full-parity-l2] indexing $(basename "${out}") ..." >&2
    samtools index "${out}"
    echo "[hc-full-parity-l2] ready ${out}" >&2
  fi
  printf '%s\n' "${out}"
}

# Only SAM paths listed in hc-full-parity case tables (not p3_malformed_* etc.).
collect_l2_sam_fixture_paths() {
  local cases_tsv rel
  for cases_tsv in "${repo_root}"/parity/fixtures/hc-full-parity/*/cases.tsv; do
    [[ -f "${cases_tsv}" ]] || continue
    while IFS= read -r rel; do
      [[ -n "${rel}" ]] && printf '%s\n' "${rel}"
    done < <(
      awk -F'\t' '
        /^#/ || NF == 0 { next }
        {
          for (i = 1; i <= NF; i++) {
            if ($i ~ /\.sam$/) print $i
          }
        }
      ' "${cases_tsv}"
    )
  done | sort -u
}

preindex_sam_fixtures_l2() {
  local rel abs
  echo "[hc-full-parity-l2] pre-indexing SAM fixtures referenced by hc-full-parity cases ..."
  while IFS= read -r rel; do
    abs="${repo_root}/${rel}"
    if [[ ! -f "${abs}" ]]; then
      echo "[hc-full-parity-l2] warn: missing SAM ${rel}" >&2
      continue
    fi
    if ! resolve_alignment_path_l2 "${abs}" >/dev/null; then
      echo "[hc-full-parity-l2] warn: skip sam->bam $(basename "${abs}")" >&2
    fi
  done < <(collect_l2_sam_fixture_paths)
  echo "[hc-full-parity-l2] SAM pre-index complete"
}

hc_dump_bin_path() {
  if [[ "${profile}" == "release" ]]; then
    printf '%s\n' "${CARGO_TARGET_DIR}/release/examples/hc_full_parity_gate_dump"
  else
    printf '%s\n' "${CARGO_TARGET_DIR}/debug/examples/hc_full_parity_gate_dump"
  fi
}

ensure_hc_dump_bin() {
  HC_DUMP_BIN="$(hc_dump_bin_path)"
  if [[ -f "${HC_DUMP_BIN}" && "${PARITY_L2_SKIP_CARGO_BUILD:-0}" == "1" ]]; then
    echo "[hc-full-parity-l2] using existing ${HC_DUMP_BIN}"
    return 0
  fi
  echo "[hc-full-parity-l2] cargo build hc_full_parity_gate_dump (profile=${profile}) ..."
  # `always` requires `width` in newer cargo; use `auto` for CI/non-TTY.
  export CARGO_TERM_PROGRESS_WHEN="${CARGO_TERM_PROGRESS_WHEN:-auto}"
  if [[ "${profile}" == "release" ]]; then
    cargo build --release -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump
  else
    cargo build -p gatk-haplotypecaller --features parity_harness --example hc_full_parity_gate_dump
  fi
  if [[ ! -x "${HC_DUMP_BIN}" ]]; then
    echo "[hc-full-parity-l2] missing binary after build: ${HC_DUMP_BIN}" >&2
    exit 127
  fi
  echo "[hc-full-parity-l2] rust dump binary ready: ${HC_DUMP_BIN}"
}
# Smoothed prob at activeProbThreshold can differ by ~2e-5 vs Java LIBS path (repeat fixture).
export PARITY_HC_ACTIVE_THRESHOLD_FUDGE="${PARITY_HC_ACTIVE_THRESHOLD_FUDGE:-2e-5}"

if [[ "${PARITY_SKIP_HC_FULL_L2:-0}" == "1" ]]; then
  echo "[hc-full-parity-l2] skipped (PARITY_SKIP_HC_FULL_L2=1)"
  exit 0
fi

profile="${PARITY_RUST_PROFILE:-dev}"

if [[ ! -f "${java_root}/b1/chr1_5_15_default_${pin_short}.tsv" ]]; then
  if [[ "${PARITY_ALLOW_MISSING_JAVA:-0}" == "1" ]]; then
    echo "[hc-full-parity-l2] skipped: no java_dumps (PARITY_ALLOW_MISSING_JAVA=1)"
    exit 0
  fi
  echo "[hc-full-parity-l2] missing java dumps; run ./scripts/parity/run_hc_full_parity_java_refresh.sh" >&2
  exit 2
fi

# Durable evidence pointers (same pattern as P12 L3/L4/L5 sign-off scripts):
# - every full-suite run tees a timestamped log
# - last_run.log always tracks the latest *full-suite* attempt (never a PARITY_L2_ONLY_PHASE slice)
# - hc_full_parity_l2_canonical.log updates only on full-suite strict green (failed=0)
l2_stamp="$(date -u '+%Y%m%dT%H%M%SZ')"
l2_log="${report_dir}/hc_full_parity_l2_${l2_stamp}.log"
l2_last_run="${report_dir}/last_run.log"
l2_canonical="${report_dir}/hc_full_parity_l2_canonical.log"
l2_full_suite=1
if [[ -n "${PARITY_L2_ONLY_PHASE:-}" ]]; then
  l2_full_suite=0
fi
if [[ "${l2_full_suite}" -eq 1 ]]; then
  exec > >(tee -a "${l2_log}") 2>&1
fi

echo "[hc-full-parity-l2] start $(date -u '+%Y-%m-%dT%H:%M:%SZ') profile=${profile} pin=${pin_short}"
if [[ "${l2_full_suite}" -eq 1 ]]; then
  echo "[hc-full-parity-l2] log ${l2_log}"
else
  echo "[hc-full-parity-l2] PARITY_L2_ONLY_PHASE=${PARITY_L2_ONLY_PHASE} (no last_run/canonical update)"
fi

echo "[hc-full-parity-l2] javac parity Java stubs ..."
if ! "${script_dir}/run_hc_full_parity_java_compile.sh"; then
  if [[ "${PARITY_ALLOW_MISSING_JAVA:-0}" == "1" ]]; then
    echo "[hc-full-parity-l2] skipped: java compile unavailable"
    exit 0
  fi
  echo "[hc-full-parity-l2] java compile failed" >&2
  exit 127
fi

ensure_hc_dump_bin
if [[ "${PARITY_L2_SKIP_SAM_PREINDEX:-0}" != "1" ]]; then
  preindex_sam_fixtures_l2
fi

passed=0
failed=0
skipped=0
l2_case_num=0
summary_json="${report_dir}/l2_summary.json"
echo "[" >"${summary_json}"
first_entry=1

record_case() {
  local phase="$1"
  local case_id="$2"
  local equal="$3"
  local rust_path="$4"
  local java_path="$5"
  local note="${6:-}"
  if [[ "${first_entry}" -eq 0 ]]; then
    echo "," >>"${summary_json}"
  fi
  first_entry=0
  PHASE="${phase}" CASE_ID="${case_id}" EQUAL="${equal}" \
    RUST_PATH="${rust_path}" JAVA_PATH="${java_path}" NOTE="${note}" \
    python3 - <<'PY' >>"${summary_json}"
import json, os
eq = os.environ["EQUAL"]
if eq == "null":
    equal = None
elif eq == "true":
    equal = True
else:
    equal = False
print(json.dumps({
    "phase": os.environ["PHASE"],
    "case_id": os.environ["CASE_ID"],
    "equal": equal,
    "rust": os.environ["RUST_PATH"],
    "java": os.environ["JAVA_PATH"],
    "note": os.environ["NOTE"],
}))
PY
}

finish_summary() {
  echo "]" >>"${summary_json}"
}

l2_phase_allowed() {
  local phase="$1"
  [[ -z "${PARITY_L2_ONLY_PHASE:-}" || "${phase}" == "${PARITY_L2_ONLY_PHASE}" ]]
}

run_l2() {
  local phase="$1"
  local case_id="$2"
  local rust_subcmd="$3"
  shift 3
  if ! l2_phase_allowed "${phase}"; then
    return 0
  fi
  local java_rel="${phase}/${case_id}_${pin_short}.tsv"
  local java_path="${java_root}/${java_rel}"
  if [[ ! -f "${java_path}" ]]; then
    echo "[hc-full-parity-l2] skip ${phase}/${case_id}: no ${java_rel}"
    skipped=$((skipped + 1))
    record_case "${phase}" "${case_id}" "null" "" "${java_path}" "missing_java_dump"
    return
  fi
  l2_case_num=$((l2_case_num + 1))
  local rust_out="${tmp_dir}/${phase}_${case_id}.rust.tsv"
  local case_json="${report_dir}/${phase}_${case_id}.json"
  local rust_stderr="${tmp_dir}/${phase}_${case_id}.rust.stderr"
  echo "[hc-full-parity-l2] [${l2_case_num}] ${phase}/${case_id} (${rust_subcmd}) ..."
  set +e
  if [[ "${PARITY_L2_VERBOSE:-0}" == "1" ]]; then
    "${HC_DUMP_BIN}" "${rust_subcmd}" "$@" >"${rust_out}" 2> >(tee "${rust_stderr}" >&2)
  else
    "${HC_DUMP_BIN}" "${rust_subcmd}" "$@" >"${rust_out}" 2>"${rust_stderr}"
  fi
  local rust_ec=$?
  set -e
  if [[ "${rust_ec}" -ne 0 ]]; then
    echo "[hc-full-parity-l2] rust dump failed ${phase}/${case_id}" >&2
    cat "${tmp_dir}/${phase}_${case_id}.rust.stderr" >&2
    failed=$((failed + 1))
    record_case "${phase}" "${case_id}" "false" "${rust_out}" "${java_path}" "rust_dump_failed"
    return
  fi
  set +e
  python3 "${compare_py}" "${rust_out}" "${java_path}" \
    --float-eps "${PARITY_L2_FLOAT_EPS:-1e-5}" \
    --json-out "${case_json}"
  # PARITY_L2_FLOAT_REL_TOL defaults to 1e-2 in compare_hc_full_parity_l2.py (LIBS vs Rust GL path).
  local cmp_ec=$?
  set -e
  if [[ "${cmp_ec}" -eq 0 ]]; then
    echo "[hc-full-parity-l2] [${l2_case_num}] ok ${phase}/${case_id}"
    passed=$((passed + 1))
    record_case "${phase}" "${case_id}" "true" "${rust_out}" "${java_path}" ""
  else
    echo "[hc-full-parity-l2] [${l2_case_num}] diff ${phase}/${case_id} (see ${case_json})"
    failed=$((failed + 1))
    record_case "${phase}" "${case_id}" "false" "${rust_out}" "${java_path}" "mismatch"
  fi
}

# B.1
while IFS=$'\t' read -r case_id ref interval padding _expected; do
  [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
  if [[ -z "${padding}" || "${padding}" == "-" ]]; then
    run_l2 b1 "${case_id}" read-shards "${repo_root}/${ref}" "${interval}"
  else
    run_l2 b1 "${case_id}" read-shards "${repo_root}/${ref}" "${interval}" "${padding}"
  fi
done <"${repo_root}/parity/fixtures/hc-full-parity/b1/cases.tsv"

run_bam_phase() {
  local phase="$1"
  local rust_subcmd="$2"
  local cases="${repo_root}/parity/fixtures/hc-full-parity/${phase}/cases.tsv"
  while IFS=$'\t' read -r case_id ref bam interval padding _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    if [[ -z "${padding}" || "${padding}" == "-" ]]; then
      run_l2 "${phase}" "${case_id}" "${rust_subcmd}" \
        "${repo_root}/${ref}" "${repo_root}/${bam}" "${interval}"
    else
      run_l2 "${phase}" "${case_id}" "${rust_subcmd}" \
        "${repo_root}/${ref}" "${repo_root}/${bam}" "${interval}" "${padding}"
    fi
  done <"${cases}"
}

run_bam_phase b2 assembly-regions

cases="${repo_root}/parity/fixtures/hc-full-parity/b2-locus/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam interval padding _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    if [[ -z "${padding}" || "${padding}" == "-" ]]; then
      run_l2 b2-locus "${case_id}" locus-pileup \
        "${repo_root}/${ref}" "${repo_root}/${bam}" "${interval}"
    else
      run_l2 b2-locus "${case_id}" locus-pileup \
        "${repo_root}/${ref}" "${repo_root}/${bam}" "${interval}" "${padding}"
    fi
  done <"${cases}"
fi
cases="${repo_root}/parity/fixtures/hc-full-parity/b2-empty-locus/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam interval padding _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    if [[ -z "${padding}" || "${padding}" == "-" ]]; then
      run_l2 b2-empty-locus "${case_id}" locus-pileup \
        "${repo_root}/${ref}" "${repo_root}/${bam}" "${interval}"
    else
      run_l2 b2-empty-locus "${case_id}" locus-pileup \
        "${repo_root}/${ref}" "${repo_root}/${bam}" "${interval}" "${padding}"
    fi
  done <"${cases}"
fi
cases="${repo_root}/parity/fixtures/hc-full-parity/b5-ref/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam interval padding _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    if [[ -z "${padding}" || "${padding}" == "-" ]]; then
      run_l2 b5-ref "${case_id}" assembly-region-reference \
        "${repo_root}/${ref}" "${repo_root}/${bam}" "${interval}"
    else
      run_l2 b5-ref "${case_id}" assembly-region-reference \
        "${repo_root}/${ref}" "${repo_root}/${bam}" "${interval}" "${padding}"
    fi
  done <"${cases}"
fi

run_bam_phase_features() {
  local cases="${repo_root}/parity/fixtures/hc-full-parity/b5-feature/cases.tsv"
  [[ -f "${cases}" ]] || return 0
  while IFS=$'\t' read -r case_id ref bam interval padding features_vcf _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    local args=("${repo_root}/${ref}" "${repo_root}/${bam}" "${interval}")
    if [[ -n "${padding}" && "${padding}" != "-" ]]; then
      args+=("${padding}")
    else
      args+=("-")
    fi
    if [[ -n "${features_vcf}" && "${features_vcf}" != "-" ]]; then
      args+=("${repo_root}/${features_vcf}")
    else
      args+=("-")
    fi
    run_l2 b5-feature "${case_id}" assembly-region-features "${args[@]}"
  done <"${cases}"
}

run_bam_phase_features

run_bam_phase_trim() {
  local cases="${repo_root}/parity/fixtures/hc-full-parity/b5-trim/cases.tsv"
  [[ -f "${cases}" ]] || return 0
  while IFS=$'\t' read -r case_id ref contig start end ext_start ext_end variants legacy _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    local args=("${repo_root}/${ref}" "${contig}" "${start}" "${end}" "${ext_start}" "${ext_end}")
    if [[ -n "${variants}" && "${variants}" != "-" ]]; then
      args+=("${repo_root}/${variants}")
    else
      args+=("-")
    fi
    if [[ "${legacy}" == "1" ]]; then
      args+=("legacy")
    fi
    run_l2 b5-trim "${case_id}" assembly-region-trim "${args[@]}"
  done <"${cases}"
}

run_bam_phase_trim

run_bam_phase_pileup_track() {
  local cases="${repo_root}/parity/fixtures/hc-full-parity/b5-pileup-track/cases.tsv"
  [[ -f "${cases}" ]] || return 0
  while IFS=$'\t' read -r case_id ref bam interval padding track _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    local args=("${repo_root}/${ref}" "${repo_root}/${bam}" "${interval}")
    if [[ -n "${padding}" && "${padding}" != "-" ]]; then
      args+=("${padding}")
    else
      args+=("-")
    fi
    if [[ "${track}" == "1" ]]; then
      args+=("1")
    fi
    run_l2 b5-pileup-track "${case_id}" assembly-region-pileup-track "${args[@]}"
  done <"${cases}"
}

run_bam_phase_pileup_track

cases="${repo_root}/parity/fixtures/hc-full-parity/b5-reads/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam interval padding _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    if [[ -z "${padding}" || "${padding}" == "-" ]]; then
      run_l2 b5-reads "${case_id}" assembly-region-reads \
        "${repo_root}/${ref}" "${repo_root}/${bam}" "${interval}"
    else
      run_l2 b5-reads "${case_id}" assembly-region-reads \
        "${repo_root}/${ref}" "${repo_root}/${bam}" "${interval}" "${padding}"
    fi
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/b5-force/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam interval padding _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    if [[ -z "${padding}" || "${padding}" == "-" ]]; then
      run_l2 b5-force "${case_id}" assembly-regions-force-active \
        "${repo_root}/${ref}" "${repo_root}/${bam}" "${interval}"
    else
      run_l2 b5-force "${case_id}" assembly-regions-force-active \
        "${repo_root}/${ref}" "${repo_root}/${bam}" "${interval}" "${padding}"
    fi
  done <"${cases}"
fi

run_bam_phase b3 apply-summary
run_bam_phase b4 walker-traversal-summary

for phase in c1 c2 c3; do
  subcmd=""
  case "${phase}" in
    c1) subcmd="raw-activity" ;;
    c2) subcmd="smoothed-activity" ;;
    c3) subcmd="active-locus" ;;
  esac
  cases="${repo_root}/parity/fixtures/hc-full-parity/${phase}/cases.tsv"
  while IFS=$'\t' read -r case_id ref bam interval _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 "${phase}" "${case_id}" "${subcmd}" \
      "${repo_root}/${ref}" "${repo_root}/${bam}" "${interval}"
  done <"${cases}"
done

cases="${repo_root}/parity/fixtures/hc-full-parity/c4-gl/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam interval _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 c4-gl "${case_id}" genotype-likelihoods \
      "${repo_root}/${ref}" "${repo_root}/${bam}" "${interval}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/c5-multi/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam interval _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 c5-multi "${case_id}" raw-activity \
      "${repo_root}/${ref}" "${repo_root}/${bam}" "${interval}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/c5-force/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam interval vcf _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 c5-force "${case_id}" raw-activity-force \
      "${repo_root}/${ref}" "${repo_root}/${bam}" "${interval}" "${repo_root}/${vcf}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/c5-ploidy/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id sample_ploidy _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 c5-ploidy "${case_id}" ploidy-resolution "${sample_ploidy}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/g1/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id pairhmm_cases _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 g1 "${case_id}" genotyping-aggregate "${repo_root}/${pairhmm_cases}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/g2-pl/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id fixture _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 g2-pl "${case_id}" genotype-format "${repo_root}/${fixture}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/g2-region/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam interval padding target _active _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
  bam_resolved="$(resolve_alignment_path_l2 "${repo_root}/${bam}")"
    args=("${repo_root}/${ref}" "${bam_resolved}" "${interval}")
    if [[ -n "${padding}" && "${padding}" != "-" ]]; then
      args+=("${padding}")
    fi
    if [[ -n "${target}" && "${target}" != "-" ]]; then
      args+=("${target}")
    fi
    run_l2 g2-region "${case_id}" assembly-region-genotype "${args[@]}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/h1/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam interval padding _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 h1 "${case_id}" reference-confidence-locus \
      "${repo_root}/${ref}" "${repo_root}/${bam}" "${interval}" "${padding}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/h1-inactive/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam interval padding _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    if [[ ! -f "${repo_root}/${ref}" ]]; then
      echo "[hc-full-parity-l2] skip h1-inactive/${case_id}: missing ref ${ref}"
      skipped=$((skipped + 1))
      continue
    fi
    bam_resolved="$(resolve_alignment_path_l2 "${repo_root}/${bam}")"
    if [[ ! -f "${bam_resolved}" ]]; then
      echo "[hc-full-parity-l2] skip h1-inactive/${case_id}: missing bam ${bam}"
      skipped=$((skipped + 1))
      continue
    fi
    run_l2 h1-inactive "${case_id}" inactive-reference-model \
      "${repo_root}/${ref}" "${bam_resolved}" "${interval}" "${padding}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/h2/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id contig length _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 h2 "${case_id}" gvcf-header "${contig}" "${length}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/h2-blocks/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id fixture _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 h2-blocks "${case_id}" gvcf-writer-blocks "${repo_root}/${fixture}"
  done <"${cases}"
fi

if [[ -f "${repo_root}/parity/fixtures/hc-full-parity/i1/expected/manifest.tsv" ]]; then
  run_l2 i1-manifest manifest annotation-manifest
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/i1/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id alt_count samples _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 i1-core "${case_id}" annotate-core "${alt_count}" "${repo_root}/${samples}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/j2/cases_vcf.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id kind ref bam_or_contig interval_or_pos ref_allele alt_allele gl_csv ad_csv _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    if [[ "${kind}" == "call-region" ]]; then
      bam_resolved="$(resolve_alignment_path_l2 "${repo_root}/${bam_or_contig}")"
      run_l2 j2-vcf "${case_id}" call-region-vcf \
        "${repo_root}/${ref}" "${bam_resolved}" "${interval_or_pos}"
    else
      run_l2 j2-vcf "${case_id}" variant-vcf-from-gl-ad \
        "${bam_or_contig}" "${interval_or_pos}" "${ref_allele}" "${alt_allele}" \
        "${gl_csv}" "${ad_csv}"
    fi
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/j2/cases_format.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id contig pos ref_allele alt_allele gl_csv ad_csv _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 j2-format "${case_id}" variant-format-from-gl-ad \
      "${contig}" "${pos}" "${ref_allele}" "${alt_allele}" "${gl_csv}" "${ad_csv}"
  done <"${cases}"
fi

while IFS=$'\t' read -r case_id bam _expected; do
  [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
  run_l2 d1 "${case_id}" read-filters "${repo_root}/${bam}"
done <"${repo_root}/parity/fixtures/hc-full-parity/d1/cases.tsv"

cases="${repo_root}/parity/fixtures/hc-full-parity/d2/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam contig pos1 cap expected mode; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    bam_path="${repo_root}/${bam}"
    if [[ -n "${mode:-}" ]]; then
      run_l2 d2 "${case_id}" downsample-positional "${bam_path}" "${cap}" "${mode}"
    else
      run_l2 d2 "${case_id}" downsample-positional "${bam_path}" "${cap}"
    fi
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/d2c/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id a b c d e f; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    if [[ "${case_id}" == target_* ]]; then
      run_l2 d2c "${case_id}" allele-biased-target-counts "${a}" "${b}"
    else
      run_l2 d2c "${case_id}" allele-biased-evidence \
        "${repo_root}/${a}" "${repo_root}/${b}" "${c}" "${d}" "${e}"
    fi
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/d2c-contam/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam interval contam _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 d2c-contam "${case_id}" raw-activity-contam \
      "${repo_root}/${ref}" "${repo_root}/${bam}" "${interval}" "${contam}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/d3/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam interval_cli _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 d3 "${case_id}" soft-clip-mean \
      "${repo_root}/${ref}" "${repo_root}/${bam}" "${interval_cli}"
  done <"${cases}"
fi

while IFS=$'\t' read -r case_id bam _expected; do
  [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
  run_l2 d4 "${case_id}" read-shard-pipeline "${repo_root}/${bam}"
done <"${repo_root}/parity/fixtures/hc-full-parity/d4/cases.tsv"

cases="${repo_root}/parity/fixtures/hc-full-parity/pre/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id bam _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 pre "${case_id}" read-pre-softclip "${repo_root}/${bam}" 0 0
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/pre-len/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id bam _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 pre-len "${case_id}" read-pre-len "${repo_root}/${bam}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/pre-mq/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id bam _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 pre-mq "${case_id}" read-pre-mq "${repo_root}/${bam}" 20
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/pre-overlap/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id bam _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 pre-overlap "${case_id}" read-pre-overlap "${repo_root}/${bam}"
  done <"${cases}"
fi

run_reads_phase() {
  local phase="$1"
  local rust_subcmd="$2"
  local cases="${repo_root}/parity/fixtures/hc-full-parity/${phase}/cases.tsv"
  [[ -f "${cases}" ]] || return 0
  while IFS=$'\t' read -r -a cols; do
    [[ -z "${cols[0]}" || "${cols[0]}" == \#* ]] && continue
    local case_id="${cols[0]}"
    local reads="${repo_root}/${cols[1]}"
  case "${phase}" in
    e1)
      run_l2 "${phase}" "${case_id}" "${rust_subcmd}" "${reads}" "${cols[2]}" "${cols[3]}"
      ;;
    e2)
      run_l2 "${phase}" "${case_id}" "${rust_subcmd}" "${reads}" "${cols[2]}" "${cols[3]}"
      ;;
    e3)
      run_l2 "${phase}" "${case_id}" "${rust_subcmd}" \
        "${reads}" "${cols[2]}" "${cols[3]}" "${cols[4]}" "${cols[5]}"
      ;;
    e4)
      run_l2 "${phase}" "${case_id}" "${rust_subcmd}" \
        "${repo_root}/${cols[1]}" "${repo_root}/${cols[2]}" \
        "${cols[3]}" "${cols[4]}" "${cols[5]}" "${cols[6]}" "${cols[7]}"
      ;;
    e5)
      run_l2 "${phase}" "${case_id}" "${rust_subcmd}" \
        "${cols[1]}" "${repo_root}/${cols[2]}" "${cols[3]}" "${cols[4]}"
      ;;
    e6)
      run_l2 "${phase}" "${case_id}" "${rust_subcmd}" \
        "${repo_root}/${cols[1]}" "${repo_root}/${cols[2]}"
      ;;
    esac
  done <"${cases}"
}

run_reads_phase e1 assembly-graph
run_reads_phase e2 assembly-graph-multi
run_reads_phase e3 assembly-graph-summary
run_reads_phase e4 assembly-graph-dangling-summary
run_reads_phase e5 assembly-graph-non-unique-summary
run_reads_phase e6 assembly-haplotype-cigars

run_reads_phase_e7() {
  local cases="${repo_root}/parity/fixtures/hc-full-parity/e7/cases.tsv"
  [[ -f "${cases}" ]] || return 0
  while IFS=$'\t' read -r case_id ref reads kmer minq min_prune min_dang recover_heads _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 e7 "${case_id}" assembly-haplotypes \
      "${repo_root}/${ref}" \
      "${repo_root}/${reads}" \
      "${kmer}" "${minq}" "${min_prune}" "${min_dang}" "${recover_heads}"
  done <"${cases}"
}

run_reads_phase_e7

run_reads_phase_e7_kbest() {
  local cases="${repo_root}/parity/fixtures/hc-full-parity/e7-kbest/cases.tsv"
  [[ -f "${cases}" ]] || return 0
  while IFS=$'\t' read -r case_id ref reads kmer minq min_prune min_dang recover_heads max_haps _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 e7-kbest "${case_id}" assembly-kbest-paths \
      "${repo_root}/${ref}" "${repo_root}/${reads}" \
      "${kmer}" "${minq}" "${min_prune}" "${min_dang}" "${recover_heads}" "${max_haps}"
  done <"${cases}"
}

run_reads_phase_e7_cap() {
  local cases="${repo_root}/parity/fixtures/hc-full-parity/e7-cap/cases.tsv"
  [[ -f "${cases}" ]] || return 0
  while IFS=$'\t' read -r case_id ref reads kmer minq min_prune min_dang recover_heads max_haps _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 e7-cap "${case_id}" assembly-haplotypes-cap \
      "${repo_root}/${ref}" "${repo_root}/${reads}" \
      "${kmer}" "${minq}" "${min_prune}" "${min_dang}" "${recover_heads}" "${max_haps}"
  done <"${cases}"
}

run_reads_phase_e7_artificial() {
  local cases="${repo_root}/parity/fixtures/hc-full-parity/e7-artificial/cases.tsv"
  [[ -f "${cases}" ]] || return 0
  while IFS=$'\t' read -r case_id ref reads kmer minq min_prune min_dang recover_heads _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 e7-artificial "${case_id}" assembly-haplotypes-production \
      "${repo_root}/${ref}" "${repo_root}/${reads}" \
      "${kmer}" "${minq}" "${min_prune}" "${min_dang}" "${recover_heads}"
  done <"${cases}"
}

run_reads_phase_e7_kbest
run_reads_phase_e7_cap
run_reads_phase_e7_artificial

run_reads_phase_e1_rec() {
  local cases="${repo_root}/parity/fixtures/hc-full-parity/e1-rec/cases.tsv"
  [[ -f "${cases}" ]] || return 0
  while IFS=$'\t' read -r case_id reads log_odds _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 e1-rec "${case_id}" read-error-correction \
      "${repo_root}/${reads}" "${log_odds}"
  done <"${cases}"
}

run_reads_phase_e8_seqgraph() {
  local cases="${repo_root}/parity/fixtures/hc-full-parity/e8-seqgraph/cases.tsv"
  [[ -f "${cases}" ]] || return 0
  while IFS=$'\t' read -r case_id ref reads kmer minq min_prune min_dang recover_heads _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 e8-seqgraph "${case_id}" assembly-seqgraph-summary \
      "${repo_root}/${ref}" "${repo_root}/${reads}" \
      "${kmer}" "${minq}" "${min_prune}" "${min_dang}" "${recover_heads}"
  done <"${cases}"
}

run_reads_phase_e0_assemble() {
  local cases="${repo_root}/parity/fixtures/hc-full-parity/e0-assemble/cases.tsv"
  [[ -f "${cases}" ]] || return 0
  while IFS=$'\t' read -r case_id ref reads _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 e0-assemble "${case_id}" assembly-assemble \
      "${repo_root}/${ref}" "${repo_root}/${reads}"
  done <"${cases}"
}

run_reads_phase_e7_junction() {
  local cases="${repo_root}/parity/fixtures/hc-full-parity/e7-junction/cases.tsv"
  [[ -f "${cases}" ]] || return 0
  while IFS=$'\t' read -r case_id ref reads kmer minq recover_edges max_haps _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 e7-junction "${case_id}" assembly-junction-haplotypes \
      "${repo_root}/${ref}" "${repo_root}/${reads}" \
      "${kmer}" "${minq}" "${recover_edges}" "${max_haps}"
  done <"${cases}"
}

run_reads_phase_e7_edges() {
  local cases="${repo_root}/parity/fixtures/hc-full-parity/e7-edges/cases.tsv"
  [[ -f "${cases}" ]] || return 0
  while IFS=$'\t' read -r case_id ref reads kmer minq recover_edges max_haps _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 e7-edges "${case_id}" assembly-junction-haplotypes \
      "${repo_root}/${ref}" "${repo_root}/${reads}" \
      "${kmer}" "${minq}" "${recover_edges}" "${max_haps}"
  done <"${cases}"
}

run_reads_phase_e1_rec
run_reads_phase_e8_seqgraph
run_reads_phase_e0_assemble
run_reads_phase_e7_junction
run_reads_phase_e7_edges

run_reads_phase_e2e() {
  local cases="${repo_root}/parity/fixtures/hc-full-parity/e2e/cases.tsv"
  [[ -f "${cases}" ]] || return 0
  while IFS=$'\t' read -r case_id ref bam interval padding target l2_strict _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    if [[ "${l2_strict}" != "yes" ]]; then
      echo "[hc-full-parity-l2] skip e2e/${case_id}: l2_strict=${l2_strict} (Rust vs Java assembly path differs; L1 only)"
      skipped=$((skipped + 1))
      record_case "e2e" "${case_id}" "null" "" "" "l2_rust_only"
      continue
    fi
    local fa="${repo_root}/${ref}"
    if [[ -f "${fa}" && ! -f "${fa}.fai" ]] && command -v samtools >/dev/null 2>&1; then
      samtools faidx "${fa}" 2>/dev/null || true
    fi
    local bam_resolved
    bam_resolved="$(resolve_alignment_path_l2 "${repo_root}/${bam}")"
    local -a extra=()
    [[ -n "${padding}" && "${padding}" != "-" ]] && extra+=("${padding}")
    [[ -n "${target}" && "${target}" != "-" && "${target}" != "active" ]] && extra+=("${target}")
    if ((${#extra[@]} > 0)); then
      run_l2 e2e "${case_id}" assembly-region-haplotypes \
        "${repo_root}/${ref}" "${bam_resolved}" "${interval}" "${extra[@]}"
    else
      run_l2 e2e "${case_id}" assembly-region-haplotypes \
        "${repo_root}/${ref}" "${bam_resolved}" "${interval}"
    fi
  done <"${cases}"
}

run_reads_phase_f1() {
  local cases="${repo_root}/parity/fixtures/hc-full-parity/f1/cases.tsv"
  [[ -f "${cases}" ]] || return 0
  while IFS=$'\t' read -r case_id cases_path _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 f1 "${case_id}" pairhmm-likelihoods "${repo_root}/${cases_path}"
  done <"${cases}"
}

run_reads_phase_f2() {
  local cases="${repo_root}/parity/fixtures/hc-full-parity/f2/cases.tsv"
  [[ -f "${cases}" ]] || return 0
  while IFS=$'\t' read -r case_id cases_path _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 f2 "${case_id}" pairhmm-native-likelihoods "${repo_root}/${cases_path}"
  done <"${cases}"
}

run_reads_phase_e2e_int() {
  local cases="${repo_root}/parity/fixtures/hc-full-parity/e2e-int/cases.tsv"
  [[ -f "${cases}" ]] || return 0
  while IFS=$'\t' read -r case_id ref bam interval padding target _l2 _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    local fa="${repo_root}/${ref}"
    if [[ -f "${fa}" && ! -f "${fa}.fai" ]] && command -v samtools >/dev/null 2>&1; then
      samtools faidx "${fa}" 2>/dev/null || true
    fi
    local bam_resolved
    bam_resolved="$(resolve_alignment_path_l2 "${repo_root}/${bam}")"
    local -a extra=()
    [[ -n "${padding}" && "${padding}" != "-" ]] && extra+=("${padding}")
    [[ -n "${target}" && "${target}" != "-" && "${target}" != "active" ]] && extra+=("${target}")
    if ((${#extra[@]} > 0)); then
      run_l2 e2e-int "${case_id}" assembly-region-pairhmm-likelihoods \
        "${repo_root}/${ref}" "${bam_resolved}" "${interval}" "${extra[@]}"
    else
      run_l2 e2e-int "${case_id}" assembly-region-pairhmm-likelihoods \
        "${repo_root}/${ref}" "${bam_resolved}" "${interval}"
    fi
  done <"${cases}"
}

run_reads_phase_f3() {
  local cases="${repo_root}/parity/fixtures/hc-full-parity/f3/cases.tsv"
  [[ -f "${cases}" ]] || return 0
  while IFS=$'\t' read -r subphase case_id cases_path _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* || "${subphase}" == \#* ]] && continue
    case "${subphase}" in
      f3-bq)
        run_l2 f3 "${case_id}" pairhmm-bq-cap "${repo_root}/${cases_path}"
        ;;
      f3-filter)
        run_l2 f3 "${case_id}" pairhmm-haplotype-filter "${repo_root}/${cases_path}"
        ;;
      *)
        echo "[hc-full-parity-l2] skip f3/${case_id}: unknown subphase ${subphase}" >&2
        skipped=$((skipped + 1))
        ;;
    esac
  done <"${cases}"
}

run_reads_phase_e2e
run_reads_phase_e2e_int
run_reads_phase_f1
run_reads_phase_f2
run_reads_phase_f3

run_reads_phase_e0_config() {
  local cases="${repo_root}/parity/fixtures/hc-full-parity/e0-config/cases.tsv"
  [[ -f "${cases}" ]] || return 0
  while IFS=$'\t' read -r case_id _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 e0-config "${case_id}" assembler-args
  done <"${cases}"
}

run_reads_phase_e0_abort() {
  local cases="${repo_root}/parity/fixtures/hc-full-parity/e0-abort/cases.tsv"
  [[ -f "${cases}" ]] || return 0
  while IFS=$'\t' read -r case_id ref reads kmer _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 e0-abort "${case_id}" assembly-graph-low-quality \
      "${repo_root}/${ref}" "${repo_root}/${reads}" "${kmer}"
  done <"${cases}"
}

run_reads_phase_d4_dragen() {
  local cases="${repo_root}/parity/fixtures/hc-full-parity/d4-dragen/cases.tsv"
  [[ -f "${cases}" ]] || return 0
  while IFS=$'\t' read -r case_id bam _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    local bam_resolved
    bam_resolved="$(resolve_alignment_path_l2 "${repo_root}/${bam}")"
    run_l2 d4-dragen "${case_id}" read-shard-pipeline-dragen "${bam_resolved}"
  done <"${cases}"
}

run_reads_phase_f4() {
  local cases="${repo_root}/parity/fixtures/hc-full-parity/f4/cases.tsv"
  [[ -f "${cases}" ]] || return 0
  while IFS=$'\t' read -r subphase case_id cases_path _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* || "${subphase}" == \#* ]] && continue
    case "${subphase}" in
      f4-config)
        run_l2 f4 "${case_id}" likelihood-engine-config
        ;;
      f4-pcr-table)
        run_l2 f4 "${case_id}" pcr-error-model
        ;;
      f4-pcr)
        run_l2 f4 "${case_id}" likelihood-pcr-read "${repo_root}/${cases_path}"
        ;;
      *)
        echo "[hc-full-parity-l2] skip f4/${case_id}: unknown subphase ${subphase}" >&2
        skipped=$((skipped + 1))
        ;;
    esac
  done <"${cases}"
}

run_reads_phase_e0_config
run_reads_phase_e0_abort
run_reads_phase_d4_dragen
run_reads_phase_f4

# --- Deferred gates L2 ---
g2af_tmp="${tmp_dir}/g2-af-inputs"
mkdir -p "${g2af_tmp}"
cases="${repo_root}/parity/fixtures/hc-full-parity/g2-af/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id gl_csv _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    input="${g2af_tmp}/${case_id}.input.tsv"
    printf '# case_id\tgl\n%s\t%s\n' "${case_id}" "${gl_csv}" >"${input}"
    run_l2 g2-af "${case_id}" af-em "${input}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/g3/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ploidy max_gt _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 g3 "${case_id}" genotype-limits "${ploidy}" "${max_gt}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/g4/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id alleles phasing phase_set _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 g4 "${case_id}" genotype-phasing "${alleles}" "${phasing}" "${phase_set}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/g4-force/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id vcf contig pos filtered _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 g4-force "${case_id}" force-calling-genotype \
      "${repo_root}/${vcf}" "${contig}" "${pos}" "${filtered}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/g-subset/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id sums is_ref max_alleles _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 g-subset "${case_id}" allele-subsetting "${sums}" "${is_ref}" "${max_alleles}"
  done <"${cases}"
fi

gsubpl_tmp="${tmp_dir}/g-subset-pl-inputs"
mkdir -p "${gsubpl_tmp}"
cases="${repo_root}/parity/fixtures/hc-full-parity/g-subset-pl/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id pl_csv ad_csv _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    input="${gsubpl_tmp}/${case_id}.input.tsv"
    printf '%s\t%s\t%s\n' "${case_id}" "${pl_csv}" "${ad_csv}" >"${input}"
    run_l2 g-subset-pl "${case_id}" subset-alleles-pl "${input}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/g-subset-vc/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id _pl ad sac _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 g-subset-vc "${case_id}" subset-alleles-vc \
      "${repo_root}/parity/fixtures/hc-full-parity/g-subset-vc/het_ac_sac_vc.tsv"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/g-subset-integration/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id sums is_ref max_alleles vc_fixture _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 g-subset-integration "${case_id}" subset-alleles-integration \
      "${sums}" "${is_ref}" "${max_alleles}" "${repo_root}/${vc_fixture}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/h2-merge/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id merge_case _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 h2-merge "${case_id}" gvcf-merge-ref-confidence "${merge_case}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/g2-subset-live/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam interval padding target max_alleles assembly_profile _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    bam_resolved="$(resolve_alignment_path_l2 "${repo_root}/${bam}")"
    args=("${repo_root}/${ref}" "${bam_resolved}" "${interval}")
    if [[ -n "${padding}" && "${padding}" != "-" ]]; then
      args+=("${padding}")
    fi
    if [[ -n "${target}" && "${target}" != "-" ]]; then
      args+=("${target}")
    fi
    args+=("${assembly_profile:--}")
    args+=("${max_alleles}")
    run_l2 g2-subset-live "${case_id}" assembly-region-genotype-subset "${args[@]}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/h2-l5/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id fixture _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 h2-l5 "${case_id}" gvcf-l5-merged "${repo_root}/${fixture}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/i1-standard/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref_fw ref_rv alt_fw alt_rv qual dp ref_bqs alt_bqs ref_pos alt_pos ref_mq alt_mq _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 i1-standard "${case_id}" standard-annotations \
      "${ref_fw}" "${ref_rv}" "${alt_fw}" "${alt_rv}" "${qual}" "${dp}" \
      "${ref_bqs}" "${alt_bqs}" "${ref_pos}" "${alt_pos}" "${ref_mq}" "${alt_mq}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/i1-as/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id site_af site_qual _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 i1-as "${case_id}" as-annotations "${site_af}" "${site_qual}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/i1-excess-het/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref het hom _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 i1-excess-het "${case_id}" excess-het "${ref}" "${het}" "${hom}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/i1-depth-hc/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ad_csv _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 i1-depth-hc "${case_id}" depth-per-sample-hc "${ad_csv}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/i1-plugins/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id plugin ref_fw ref_rv alt_fw alt_rv qual dp ref_bqs alt_bqs _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 i1-plugins "${case_id}" annotation-plugin \
      "${plugin}" "${ref_fw}" "${ref_rv}" "${alt_fw}" "${alt_rv}" \
      "${qual}" "${dp}" "${ref_bqs}" "${alt_bqs}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/j-modes/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id mode has_variant locus_count _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 j-modes "${case_id}" emit-mode-decision "${mode}" "${has_variant}" "${locus_count}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/j-bamout/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id enabled write_count _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 j-bamout "${case_id}" bamout-stub "${enabled}" "${write_count}"
  done <"${cases}"
fi

run_l2 j-dragen default_off dragen-mode-branch

cases="${repo_root}/parity/fixtures/hc-full-parity/pre-dragstr/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id params_loaded _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 pre-dragstr "${case_id}" dragstr-calibration "${params_loaded}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/e-debug/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id failure_bam graph_dot _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    run_l2 e-debug "${case_id}" assembly-debug-stub "${failure_bam}" "${graph_dot}"
  done <"${cases}"
fi

finish_summary

echo "[hc-full-parity-l2] done $(date -u '+%Y-%m-%dT%H:%M:%SZ') passed=${passed} failed=${failed} skipped=${skipped} cases_run=${l2_case_num} strict=${strict}"
echo "[hc-full-parity-l2] summary ${summary_json}"

if [[ "${l2_full_suite}" -eq 1 ]]; then
  # `exec > >(tee …)` is asynchronous — wait until the done line is on disk before
  # copying evidence pointers so last_run/canonical are not truncated mid-flush.
  for _ in $(seq 1 100); do
    if grep -qF "passed=${passed} failed=${failed}" "${l2_log}" 2>/dev/null; then
      break
    fi
    sleep 0.05
  done
  # Always point last_run at the newest full-suite attempt (green or red).
  cp -f "${l2_log}" "${l2_last_run}"
  if [[ "${failed}" -eq 0 && "${strict}" == "1" && "${l2_case_num}" -gt 0 ]]; then
    cp -f "${l2_log}" "${l2_canonical}"
    echo "[hc-full-parity-l2] last_run ${l2_last_run}"
    echo "[hc-full-parity-l2] canonical ${l2_canonical} (full-suite strict green)"
  else
    echo "[hc-full-parity-l2] last_run ${l2_last_run}"
    echo "[hc-full-parity-l2] canonical unchanged (failed=${failed} strict=${strict} cases_run=${l2_case_num})"
  fi
fi

if [[ "${failed}" -gt 0 && "${strict}" == "1" ]]; then
  exit 1
fi
# Advisory default: mismatches expected until gap closure (see JAVA_RUST_ALGORITHM_GAPS.md).
exit 0
