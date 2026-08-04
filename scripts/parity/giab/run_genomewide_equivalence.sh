#!/usr/bin/env bash
# GIAB multi-sample HaplotypeCaller equivalence: Java GATK4 (pinned) vs gatk-rs,
# evaluated with gatk-rs-equiv (hap.py / RTG) + /usr/bin/time resource capture.
#
# Modes (GIAB_MODE) — what “genome-wide” means here:
#   smoke       small windows only (M4 / PR)
#   ci-subset   FULL chr20+chr21 + 50kb probes on other autosomes  ← CI default
#   chr20-21    full chr20+chr21 only
#   autosomes   full chr1–22 (big iron only)
#
# Phases (GIAB_PHASE) — for GitHub-hosted 6h job limit:
#   all       prepare + HC shards + concat + equiv + dashboard (default; local)
#   prepare   intervals/shards/SCOPE only (no HC)
#   hc        run filtered HC shard×engine jobs (no concat/equiv)
#   finalize  concat shard VCFs + equiv + dashboard (no HC)
#
# Filters (hc phase):
#   GIAB_HC_SHARDS=00_chr20_w00      comma-separated shard basenames (empty = all)
#   GIAB_HC_ENGINES=java|rust|java,rust
#   GIAB_HC_WINDOW_BP=1000000        window size for splitting full chr20/21 shards
#
# Usage:
#   GIAB_MODE=smoke GIAB_SAMPLES=HG001 ./scripts/parity/giab/run_genomewide_equivalence.sh
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../../.." && pwd)"
# shellcheck source=../m4_disk_guard.sh
source "${repo_root}/scripts/parity/m4_disk_guard.sh"
# shellcheck source=lib_giab.sh
source "${repo_root}/scripts/parity/giab/lib_giab.sh"

phase="${GIAB_PHASE:-all}"
case "${phase}" in
  all|prepare|hc|finalize) ;;
  *)
    echo "[giab] unknown GIAB_PHASE=${phase} (want all|prepare|hc|finalize)" >&2
    exit 2
    ;;
esac

min_free_default=12
if [[ "${phase}" == "prepare" ]]; then
  min_free_default=6
fi
m4_require_free_gb "${GIAB_RUN_MIN_FREE_GB:-${min_free_default}}" || exit 1

mode="${GIAB_MODE:-ci-subset}"
samples_csv="${GIAB_SAMPLES:-HG001}"
out_root="${GIAB_OUT_ROOT:-${repo_root}/parity/giab/runs/$(date -u +%Y%m%dT%H%M%SZ)_${mode}}"
truth_root="${GIAB_TRUTH_ROOT:-${repo_root}/parity/giab/truth}"
strat_root="${GIAB_STRAT_ROOT:-${repo_root}/parity/giab/stratifications/GRCh37}"
ref="${GIAB_REFERENCE:-${repo_root}/parity/realworld/assets/hs37d5.simple.fa}"
rust_bin="${GIAB_RUST_BIN:-${CARGO_TARGET_DIR:-${repo_root}/target}/release/gatk-rs}"
equiv_bin="${GIAB_EQUIV_BIN:-${CARGO_TARGET_DIR:-${repo_root}/target}/release/gatk-rs-equiv}"
f1_delta="${GIAB_F1_DELTA_THRESHOLD:-0.02}"
threads="${GIAB_THREADS:-2}"
fetch_truth="${GIAB_FETCH_TRUTH:-1}"
stage_ref="${GIAB_STAGE_REF:-1}"
java_jar="${GIAB_JAVA_GATK_JAR:-}"
java_bin="${GIAB_JAVA_GATK_BIN:-}"
skip_equiv_engine="${GIAB_SKIP_EQUIV_ENGINE:-0}"
hc_shards_csv="${GIAB_HC_SHARDS:-}"
hc_engines_csv="${GIAB_HC_ENGINES:-java,rust}"
ftp_data="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data"

if [[ "${phase}" == "hc" ]]; then
  skip_equiv_engine=1
fi

bam_url_for_sample() {
  case "$1" in
    HG001) echo "${GIAB_HG001_BAM_URL:-${ftp_data}/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam}" ;;
    HG002) echo "${GIAB_HG002_BAM_URL:-${ftp_data}/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.300x.bam}" ;;
    HG005) echo "${GIAB_HG005_BAM_URL:-${ftp_data}/ChineseTrio/HG005_NA24631_son/HG005_NA24631_son_HiSeq_300x/NHGRI_Illumina300X_Chinesetrio_novoalign_bams/HG005.hs37d5.300x.bam}" ;;
    *) return 1 ;;
  esac
}

csv_contains() {
  # csv_contains haystack_csv needle
  local csv="$1" needle="$2" item
  local OLD_IFS="${IFS}"
  IFS=','
  # shellcheck disable=SC2086
  set -- ${csv}
  IFS="${OLD_IFS}"
  for item in "$@"; do
    item="$(echo "${item}" | tr -d '[:space:]')"
    [[ "${item}" == "${needle}" ]] && return 0
  done
  return 1
}

want_shard() {
  local name="$1"
  [[ -z "${hc_shards_csv}" ]] && return 0
  csv_contains "${hc_shards_csv}" "${name}"
}

want_engine() {
  local eng="$1"
  csv_contains "${hc_engines_csv}" "${eng}"
}

mkdir -p "${out_root}"
echo "=== run_genomewide_equivalence ==="
echo "phase=${phase}"
echo "mode=${mode}"
echo "scope: $(giab_mode_description "${mode}")"
echo "samples=${samples_csv}"
echo "out_root=${out_root}"
echo "hc_shards=${hc_shards_csv:-<all>}"
echo "hc_engines=${hc_engines_csv}"

if [[ "${phase}" == "all" || "${phase}" == "prepare" ]]; then
  giab_mode_description "${mode}" > "${out_root}/SCOPE.txt"
  {
    echo "GIAB_MODE=${mode}"
    echo "GIAB_SAMPLES=${samples_csv}"
    echo "created_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "description=$(giab_mode_description "${mode}")"
  } > "${out_root}/run.env"
fi

if [[ "${fetch_truth}" == "1" ]]; then
  "${repo_root}/scripts/parity/giab/fetch_giab_truthsets.sh"
fi

if [[ "${stage_ref}" == "1" && ! -f "${ref}" ]]; then
  echo "[giab] staging reference via realworld step 03…"
  "${repo_root}/scripts/parity/realworld/03_stage_reference_and_truth.sh"
  ref="${repo_root}/parity/realworld/assets/hs37d5.simple.fa"
fi
if [[ ! -f "${ref}" ]]; then
  echo "[giab] missing reference: ${ref}" >&2
  if [[ "${stage_ref}" != "1" ]]; then
    echo "[giab] GIAB_STAGE_REF=${stage_ref}: expected hs37d5 from prepare artifact / cache; refusing to curl EBI FTP." >&2
  fi
  exit 2
fi
if [[ ! -f "${ref}.fai" ]]; then
  echo "[giab] missing reference index: ${ref}.fai" >&2
  exit 2
fi
ref_dict="$(dirname "${ref}")/hs37d5.simple.dict"
if [[ ! -f "${ref_dict}" ]]; then
  echo "[giab] missing reference dict: ${ref_dict}" >&2
  exit 2
fi

if [[ "${phase}" != "prepare" ]]; then
  if [[ ! -x "${rust_bin}" ]]; then
    echo "[giab] building gatk-rs release…"
    "${repo_root}/scripts/parity/build_gatk_rs_release.sh"
    rust_bin="${CARGO_TARGET_DIR:-${repo_root}/target}/release/gatk-rs"
  fi
  if [[ ! -x "${equiv_bin}" && "${skip_equiv_engine}" != "1" ]]; then
    echo "[giab] building gatk-rs-equiv release…"
    cargo build -p gatk-rs-equiv --release -j "${CARGO_BUILD_JOBS:-1}"
    equiv_bin="${CARGO_TARGET_DIR:-${repo_root}/target}/release/gatk-rs-equiv"
  fi
fi

if [[ "${phase}" != "prepare" ]] && ! command -v samtools >/dev/null 2>&1; then
  echo "[giab] samtools required on PATH" >&2
  exit 2
fi

if [[ "${phase}" == "all" || "${phase}" == "prepare" || ! -f "${out_root}/intervals.txt" ]]; then
  giab_build_intervals "${mode}" > "${out_root}/intervals.txt"
fi
echo "[giab] intervals:"
cat "${out_root}/intervals.txt"

# Load full-mode intervals (finalize BAM staging / smoke).
intervals=()
while IFS= read -r line || [[ -n "${line}" ]]; do
  [[ -z "${line}" ]] && continue
  intervals+=("${line}")
done < "${out_root}/intervals.txt"

shard_root="${out_root}/shards"
if [[ "${phase}" == "all" || "${phase}" == "prepare" || ! -d "${shard_root}" ]]; then
  giab_write_hc_shards "${shard_root}" "${out_root}/intervals.txt" "${mode}"
fi
echo "[giab] HC shards:"
ls -1 "${shard_root}"/*.intervals 2>/dev/null | while IFS= read -r sf; do
  echo "  $(basename "${sf}"): $(wc -l < "${sf}" | tr -d ' ') interval(s)"
done

python3 - "${shard_root}" "${out_root}/shards.json" <<'PY'
import json, pathlib, sys
shard_root = pathlib.Path(sys.argv[1])
out = pathlib.Path(sys.argv[2])
shards = sorted(p.stem for p in shard_root.glob("*.intervals"))
out.write_text(json.dumps(shards), encoding="utf-8")
print(f"[giab] wrote {out} → {shards}")
PY

if [[ "${phase}" == "prepare" ]]; then
  echo "[giab] prepare phase complete out=${out_root}"
  exit 0
fi

time_backend="$(giab_time_backend)"
echo "[giab] time_backend=${time_backend}"

overall_gate=0
summary_jsonl="${out_root}/samples.jsonl"
if [[ "${phase}" != "hc" ]]; then
  : > "${summary_jsonl}"
fi

stage_bam_for_interval_list() {
  local sample="$1" bam_url="$2" out_bam="$3"
  shift 3
  local -a ivs=("$@")
  local bai_local="${out_bam}.remote.bai"
  local out_bai="${out_bam}.bai"
  mkdir -p "$(dirname "${out_bam}")"
  if [[ -f "${out_bam}" && -f "${out_bai}" ]]; then
    echo "[giab] reuse BAM ${out_bam}"
    return 0
  fi
  if [[ "${#ivs[@]}" -eq 0 ]]; then
    echo "[giab] stage_bam: empty interval list" >&2
    return 2
  fi
  if [[ ! -f "${bai_local}" ]]; then
    echo "[giab] download BAI for ${sample}…"
    curl -fL --retry 3 -o "${bai_local}.partial" "${bam_url}.bai"
    mv -f "${bai_local}.partial" "${bai_local}"
  fi
  echo "[giab] slicing remote BAM for ${sample} (${#ivs[@]} intervals)…"
  samtools view -b -X "${bam_url}" "${bai_local}" "${ivs[@]}" > "${out_bam}.partial"
  mv -f "${out_bam}.partial" "${out_bam}"
  samtools index "${out_bam}" "${out_bai}"
}

load_intervals_file() {
  # load_intervals_file path → fills loaded_intervals array
  local path="$1" line
  loaded_intervals=()
  while IFS= read -r line || [[ -n "${line}" ]]; do
    [[ -z "${line}" ]] && continue
    loaded_intervals+=("${line}")
  done < "${path}"
}

run_timed() {
  local label="$1" log="$2"
  shift 2
  echo "[giab] RUN ${label}: $*"
  giab_run_timed "${time_backend}" "${log}" "$@"
}

OLD_IFS="${IFS}"
IFS=','
# shellcheck disable=SC2086
set -- ${samples_csv}
IFS="${OLD_IFS}"

for sample in "$@"; do
  sample="$(echo "${sample}" | tr -d '[:space:]')"
  [[ -n "${sample}" ]] || continue
  echo "======== sample ${sample} ========"
  sdir="${out_root}/${sample}"
  mkdir -p "${sdir}/hc" "${sdir}/equiv" "${sdir}/time"

  truth_vcf="${truth_root}/${sample}_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
  truth_bed="${truth_root}/${sample}_GRCh37_1_22_v4.2.1_benchmark.bed"
  if [[ ! -f "${truth_vcf}" || ! -f "${truth_bed}" ]]; then
    echo "[giab] missing truth for ${sample} under ${truth_root}" >&2
    exit 2
  fi

  bam_url="$(bam_url_for_sample "${sample}")" || {
    echo "[giab] no BAM URL for ${sample}" >&2
    exit 2
  }

  java_vcf="${sdir}/hc/java.vcf"
  rust_vcf="${sdir}/hc/rust.vcf"
  java_time="${sdir}/time/java.time.txt"
  rust_time="${sdir}/time/rust.time.txt"
  shard_vcf_dir="${sdir}/hc/shards"
  mkdir -p "${shard_vcf_dir}"

  run_hc_shard() {
    local engine="$1" shard_name="$2" shard_intervals="$3" out_vcf="$4" time_log="$5" bam_path="$6"
    local l_args=()
    local iv
    while IFS= read -r iv || [[ -n "${iv}" ]]; do
      [[ -z "${iv}" ]] && continue
      l_args+=(-L "${iv}")
    done < "${shard_intervals}"
    if [[ "${#l_args[@]}" -eq 0 ]]; then
      echo "[giab] empty shard ${shard_name}; skipping" >&2
      return 0
    fi
    if [[ "${engine}" == "java" ]]; then
      if [[ -n "${java_jar}" ]]; then
        run_timed "java-hc-${sample}-${shard_name}" "${time_log}" \
          java -Xmx4g -jar "${java_jar}" HaplotypeCaller \
            -R "${ref}" -I "${bam_path}" -O "${out_vcf}" --verbosity ERROR \
            --native-pair-hmm-threads "${threads}" "${l_args[@]}"
      elif [[ -n "${java_bin}" ]]; then
        run_timed "java-hc-${sample}-${shard_name}" "${time_log}" \
          "${java_bin}" HaplotypeCaller \
            -R "${ref}" -I "${bam_path}" -O "${out_vcf}" --verbosity ERROR \
            --native-pair-hmm-threads "${threads}" "${l_args[@]}"
      elif command -v gatk >/dev/null 2>&1; then
        run_timed "java-hc-${sample}-${shard_name}" "${time_log}" \
          gatk HaplotypeCaller \
            -R "${ref}" -I "${bam_path}" -O "${out_vcf}" --verbosity ERROR \
            --native-pair-hmm-threads "${threads}" "${l_args[@]}"
      else
        # shellcheck source=../lib_pinned_gatk.sh
        source "${repo_root}/scripts/parity/lib_pinned_gatk.sh"
        local img="${GATK_DOCKER_IMAGE:-us.gcr.io/broad-gatk/gatk:4.4.0.0}"
        local plat="${GATK_DOCKER_PLATFORM:-linux/amd64}"
        run_timed "java-hc-${sample}-${shard_name}" "${time_log}" \
          docker run --rm --platform "${plat}" \
            -v "${repo_root}:${repo_root}" -w "${repo_root}" \
            "${img}" gatk HaplotypeCaller \
            -R "${ref}" -I "${bam_path}" -O "${out_vcf}" --verbosity ERROR \
            --native-pair-hmm-threads "${threads}" "${l_args[@]}"
      fi
    else
      run_timed "rust-hc-${sample}-${shard_name}" "${time_log}" \
        env RAYON_NUM_THREADS="${threads}" \
        "${rust_bin}" HaplotypeCaller \
          -R "${ref}" -I "${bam_path}" -O "${out_vcf}" "${l_args[@]}"
    fi
  }

  if [[ "${phase}" == "all" || "${phase}" == "hc" ]]; then
    for shard_file in "${shard_root}"/*.intervals; do
      [[ -f "${shard_file}" ]] || continue
      shard_name="$(basename "${shard_file}" .intervals)"
      want_shard "${shard_name}" || continue

      load_intervals_file "${shard_file}"
      if [[ "${mode}" == "smoke" ]]; then
        # Distinct filename so older full-30× P12 smoke slices are not reused.
        bam="${sdir}/hc/${sample}.${mode}.${shard_name}.p12_20k.bam"
        giab_stage_smoke_bam_hybrid "${sample}" "${bam_url}" "${bam}" "${repo_root}" "${loaded_intervals[@]}"
      else
        bam="${sdir}/hc/${sample}.${mode}.${shard_name}.bam"
        stage_bam_for_interval_list "${sample}" "${bam_url}" "${bam}" "${loaded_intervals[@]}"
      fi

      java_shard_vcf="${shard_vcf_dir}/java.${shard_name}.vcf"
      rust_shard_vcf="${shard_vcf_dir}/rust.${shard_name}.vcf"
      java_shard_time="${sdir}/time/java.${shard_name}.time.txt"
      rust_shard_time="${sdir}/time/rust.${shard_name}.time.txt"

      if want_engine java; then
        if [[ ! -f "${java_shard_vcf}" || "${GIAB_FORCE_HC:-0}" == "1" ]]; then
          run_hc_shard java "${shard_name}" "${shard_file}" "${java_shard_vcf}" "${java_shard_time}" "${bam}"
        else
          echo "[giab] reuse ${java_shard_vcf}"
        fi
      fi
      if want_engine rust; then
        if [[ ! -f "${rust_shard_vcf}" || "${GIAB_FORCE_HC:-0}" == "1" ]]; then
          run_hc_shard rust "${shard_name}" "${shard_file}" "${rust_shard_vcf}" "${rust_shard_time}" "${bam}"
        else
          echo "[giab] reuse ${rust_shard_vcf}"
        fi
      fi
    done
  fi

  if [[ "${phase}" == "hc" ]]; then
    echo "[giab] hc phase complete for ${sample}"
    continue
  fi

  # finalize / all: concat every shard VCF present for this sample.
  java_shard_vcfs=()
  rust_shard_vcfs=()
  for shard_file in "${shard_root}"/*.intervals; do
    [[ -f "${shard_file}" ]] || continue
    shard_name="$(basename "${shard_file}" .intervals)"
    java_shard_vcf="${shard_vcf_dir}/java.${shard_name}.vcf"
    rust_shard_vcf="${shard_vcf_dir}/rust.${shard_name}.vcf"
    if [[ ! -f "${java_shard_vcf}" ]]; then
      echo "[giab] missing Java shard VCF: ${java_shard_vcf}" >&2
      exit 2
    fi
    if [[ ! -f "${rust_shard_vcf}" ]]; then
      echo "[giab] missing Rust shard VCF: ${rust_shard_vcf}" >&2
      exit 2
    fi
    java_shard_vcfs+=("${java_shard_vcf}")
    rust_shard_vcfs+=("${rust_shard_vcf}")
  done

  if [[ "${#java_shard_vcfs[@]}" -eq 0 ]]; then
    echo "[giab] no HC shards produced for ${sample}" >&2
    exit 2
  fi

  echo "[giab] concat ${#java_shard_vcfs[@]} Java shard VCF(s) → ${java_vcf}"
  giab_concat_vcfs "${java_vcf}" "${java_shard_vcfs[@]}"
  : > "${java_time}"
  echo "[giab] concat ${#rust_shard_vcfs[@]} Rust shard VCF(s) → ${rust_vcf}"
  giab_concat_vcfs "${rust_vcf}" "${rust_shard_vcfs[@]}"
  : > "${rust_time}"

  # Equiv may inspect the BAM; stage the full interval union.
  if [[ "${mode}" == "smoke" ]]; then
    bam="${sdir}/hc/${sample}.${mode}.p12_20k.bam"
    giab_stage_smoke_bam_hybrid "${sample}" "${bam_url}" "${bam}" "${repo_root}" "${intervals[@]}"
  else
    bam="${sdir}/hc/${sample}.${mode}.bam"
    stage_bam_for_interval_list "${sample}" "${bam_url}" "${bam}" "${intervals[@]}"
  fi

  java_perf="$(giab_parse_time_log "${java_time}" 2>/dev/null || echo '{}')"
  rust_perf="$(giab_parse_time_log "${rust_time}" 2>/dev/null || echo '{}')"

  equiv_rc=0
  if [[ "${skip_equiv_engine}" != "1" ]]; then
    mkdir -p "${sdir}/equiv"
    # Prefer hardlink (same inode, no duplex) when possible; fall back to cp.
    if ! ln -f "${java_vcf}" "${sdir}/equiv/java.vcf" 2>/dev/null; then
      cp -f "${java_vcf}" "${sdir}/equiv/java.vcf"
    fi
    if ! ln -f "${rust_vcf}" "${sdir}/equiv/rust.vcf" 2>/dev/null; then
      cp -f "${rust_vcf}" "${sdir}/equiv/rust.vcf"
    fi
    # After successful concat, drop per-shard VCFs unless retained for debugging.
    if [[ "${GIAB_KEEP_SHARD_VCFS:-0}" != "1" ]]; then
      rm -f "${java_shard_vcfs[@]}" "${rust_shard_vcfs[@]}" 2>/dev/null || true
    fi
    strat_args=()
    mkdir -p "${sdir}/equiv/strat"
    for bedgz in \
      "${strat_root}/GRCh37_AllTandemRepeatsandHomopolymers_slop5.bed.gz" \
      "${strat_root}/GRCh37_segdups.bed.gz" \
      "${strat_root}/GRCh37_alldifficultregions.bed.gz"
    do
      [[ -f "${bedgz}" ]] || continue
      base="$(basename "${bedgz}" .bed.gz)"
      gunzip -c "${bedgz}" > "${sdir}/equiv/strat/${base}.bed"
      name="${base#GRCh37_}"
      strat_args+=(--stratification-bed "${name}=${sdir}/equiv/strat/${base}.bed")
    done

    set +e
    "${equiv_bin}" run \
      --rust-binary "${rust_bin}" \
      --reference "${ref}" \
      --bam "${bam}" \
      --truth-vcf "${truth_vcf}" \
      --confident-regions "${truth_bed}" \
      --out "${sdir}/equiv" \
      --reuse-vcfs \
      --threads "${threads}" \
      --f1-delta-threshold "${f1_delta}" \
      --min-free-gb "${GIAB_EQUIV_MIN_FREE_GB:-6}" \
      "${strat_args[@]}"
    equiv_rc=$?
    set -e

    # Free finalize BAM + gunzipped strat beds after scoring (Air / limited-HDD recipe).
    if [[ "${GIAB_KEEP_EQUIV_INTERMEDIATES:-0}" != "1" ]]; then
      rm -f "${bam}" "${bam}.bai" 2>/dev/null || true
      rm -rf "${sdir}/equiv/strat" 2>/dev/null || true
    fi
  else
    echo "[giab] GIAB_SKIP_EQUIV_ENGINE=1 — skipping hap.py/RTG"
  fi

  if [[ "${equiv_rc}" -ne 0 ]]; then
    overall_gate=1
  fi

  python3 - "${summary_jsonl}" "${sample}" "${mode}" "${equiv_rc}" "${java_perf}" "${rust_perf}" "${sdir}" <<'PY'
import json, pathlib, sys
path, sample, mode, rc, jp, rp, sdir = sys.argv[1:8]
row = {
    "sample": sample,
    "mode": mode,
    "equiv_exit": int(rc),
    "gate_passed": int(rc) == 0,
    "java_perf": json.loads(jp or "{}"),
    "rust_perf": json.loads(rp or "{}"),
    "dir": sdir,
}
results = pathlib.Path(sdir) / "equiv" / "results.json"
if results.is_file():
    row["equiv_results"] = json.loads(results.read_text(encoding="utf-8"))
with open(path, "a", encoding="utf-8") as fh:
    fh.write(json.dumps(row) + "\n")
PY
done

if [[ "${phase}" == "hc" ]]; then
  echo "[giab] hc phase finished out=${out_root}"
  exit 0
fi

python3 "${repo_root}/scripts/parity/giab/build_dashboard.py" \
  --run-dir "${out_root}" \
  --out-dir "${out_root}/dashboard"

echo "[giab] finished overall_gate=${overall_gate} out=${out_root}"
echo "[giab] SCOPE: $(giab_mode_description "${mode}")"
exit "${overall_gate}"
