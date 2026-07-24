#!/usr/bin/env bash
# Regenerate frozen Java L2 dumps under parity/fixtures/hc-full-parity/java_dumps/.
#
# A.0.8 contract: keep this script aligned with `run_hc_full_parity_l2.sh`:
#   every phase/case that strict L2 compares must have a matching `write_dump` loop
#   here (same subcommand + paths as the Rust driver).
set -euo pipefail

script_dir="$(cd "$(dirname "$0")" && pwd)"
# shellcheck source=lib_pinned_gatk.sh
source "${script_dir}/lib_pinned_gatk.sh"

repo_root="${GATK_RS_REPO_ROOT}"
cd "${repo_root}"

pin_short="${GATK_PINNED_SHA_SHORT:-2dbc0258}"
out_root="${repo_root}/parity/fixtures/hc-full-parity/java_dumps"
dump_sh="${script_dir}/run_hc_full_parity_java_dump.sh"

if [[ "${PARITY_SKIP_HC_FULL_JAVA_REFRESH:-0}" == "1" ]]; then
  echo "[hc-full-parity-java-refresh] skipped"
  exit 0
fi

mkdir -p "${out_root}"

maybe_index_bam() {
  local bam="$1"
  [[ "${bam}" == *.bam ]] || return 0
  if [[ -f "${bam}" && ! -f "${bam}.bai" && ! -f "${bam%.bam}.bai" ]]; then
    if command -v samtools >/dev/null 2>&1; then
      echo "[hc-full-parity-java-refresh] indexing ${bam}"
      samtools index "${bam}"
    fi
  fi
}

maybe_index_fa() {
  local fa="$1"
  if [[ ! -f "${fa}" ]]; then
    return
  fi
  if [[ ! -f "${fa}.fai" ]] && command -v samtools >/dev/null 2>&1; then
    echo "[hc-full-parity-java-refresh] faidx ${fa}" >&2
    samtools faidx "${fa}"
  fi
  local dict="${fa%.fa}.dict"
  if [[ ! -f "${dict}" && ! -f "${fa%.fa}.dict" ]]; then
    echo "[hc-full-parity-java-refresh] CreateSequenceDictionary ${fa}" >&2
    "${script_dir}/run_java_gatk.sh" /dev/null CreateSequenceDictionary -R "${fa}" </dev/null || true
  fi
}

# GATK interval traversal requires indexed BAM; convert SAM fixtures on demand.
resolve_alignment_path() {
  local path="$1"
  if [[ "${path}" != *.sam ]]; then
    echo "${path}"
    return
  fi
  local cache_dir="${repo_root}/parity/build/sam-indexed-bam"
  mkdir -p "${cache_dir}"
  local out="${cache_dir}/$(basename "${path%.sam}").bam"
  if [[ ! -f "${out}" ]]; then
    echo "[hc-full-parity-java-refresh] sam->bam ${path}" >&2
    samtools view -bS "${path}" | samtools sort -o "${out}"
    samtools index "${out}"
  fi
  printf '%s\n' "${out}"
}

write_dump() {
  local rel="$1"
  shift
  local out="${out_root}/${rel}"
  mkdir -p "$(dirname "${out}")"
  echo "[hc-full-parity-java-refresh] ${rel}"
  set +e
  "${dump_sh}" "$@" >"${out}" 2>"${out}.stderr" </dev/null
  local ec=$?
  set -e
  if [[ "${ec}" -ne 0 ]]; then
    echo "[hc-full-parity-java-refresh] dump failed ec=${ec} (${rel})" >&2
    cat "${out}.stderr" >&2
    exit "${ec}"
  fi
}

# B.1
while IFS=$'\t' read -r case_id ref interval padding _expected; do
  [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
  maybe_index_fa "${repo_root}/${ref}"
  pad_arg="${padding}"
  [[ -z "${pad_arg}" || "${pad_arg}" == "-" ]] && pad_arg="-"
  write_dump "b1/${case_id}_${pin_short}.tsv" read-shards "${repo_root}/${ref}" "${interval}" "${pad_arg}"
done <"${repo_root}/parity/fixtures/hc-full-parity/b1/cases.tsv"

# B.2.1 locus pileup depth
cases="${repo_root}/parity/fixtures/hc-full-parity/b2-locus/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam interval padding _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    maybe_index_fa "${repo_root}/${ref}"
    bam_resolved="$(resolve_alignment_path "${repo_root}/${bam}" | tail -1)"
    maybe_index_bam "${bam_resolved}"
    pad="${padding:--}"
    write_dump "b2-locus/${case_id}_${pin_short}.tsv" \
      locus-pileup "${repo_root}/${ref}" "${bam_resolved}" "${interval}" "${pad}"
  done <"${cases}"
fi

# B.2.2 empty locus / zero-depth pileup (same subcommand as b2-locus; separate L2 phase for manifest)
cases="${repo_root}/parity/fixtures/hc-full-parity/b2-empty-locus/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam interval padding _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    maybe_index_fa "${repo_root}/${ref}"
    bam_resolved="$(resolve_alignment_path "${repo_root}/${bam}" | tail -1)"
    maybe_index_bam "${bam_resolved}"
    pad="${padding:--}"
    write_dump "b2-empty-locus/${case_id}_${pin_short}.tsv" \
      locus-pileup "${repo_root}/${ref}" "${bam_resolved}" "${interval}" "${pad}"
  done <"${cases}"
fi

# B.5.8 pileup tracking on regions
cases="${repo_root}/parity/fixtures/hc-full-parity/b5-pileup-track/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam interval padding track _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    maybe_index_fa "${repo_root}/${ref}"
    bam_resolved="$(resolve_alignment_path "${repo_root}/${bam}" | tail -1)"
    maybe_index_bam "${bam_resolved}"
    pad="${padding:--}"
    pileup_args=(
      assembly-region-pileup-track "${repo_root}/${ref}" "${bam_resolved}" "${interval}" "${pad}"
    )
    [[ "${track}" == "1" ]] && pileup_args+=("1")
    write_dump "b5-pileup-track/${case_id}_${pin_short}.tsv" "${pileup_args[@]}"
  done <"${cases}"
fi

# B.5.5 AssemblyRegionTrimmer
cases="${repo_root}/parity/fixtures/hc-full-parity/b5-trim/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref contig start end ext_start ext_end variants legacy _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    maybe_index_fa "${repo_root}/${ref}"
    if [[ -n "${variants}" && "${variants}" != "-" ]]; then
      v="${repo_root}/${variants}"
    else
      v="-"
    fi
    trim_args=(
      assembly-region-trim "${repo_root}/${ref}" "${contig}" "${start}" "${end}"
      "${ext_start}" "${ext_end}" "${v}"
    )
    [[ "${legacy}" == "1" ]] && trim_args+=("legacy")
    write_dump "b5-trim/${case_id}_${pin_short}.tsv" "${trim_args[@]}"
  done <"${cases}"
fi

# B.5.4 FeatureContext per region
cases="${repo_root}/parity/fixtures/hc-full-parity/b5-feature/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam interval padding features_vcf _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    maybe_index_fa "${repo_root}/${ref}"
    bam_resolved="$(resolve_alignment_path "${repo_root}/${bam}" | tail -1)"
    maybe_index_bam "${bam_resolved}"
    pad="${padding:--}"
    if [[ -n "${features_vcf}" && "${features_vcf}" != "-" ]]; then
      feat="${repo_root}/${features_vcf}"
    else
      feat="-"
    fi
    write_dump "b5-feature/${case_id}_${pin_short}.tsv" \
      assembly-region-features "${repo_root}/${ref}" "${bam_resolved}" "${interval}" "${pad}" "${feat}"
  done <"${cases}"
fi

# B.5.3 ReferenceContext per region
cases="${repo_root}/parity/fixtures/hc-full-parity/b5-ref/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam interval padding _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    maybe_index_fa "${repo_root}/${ref}"
    bam_resolved="$(resolve_alignment_path "${repo_root}/${bam}" | tail -1)"
    maybe_index_bam "${bam_resolved}"
    pad="${padding:--}"
    write_dump "b5-ref/${case_id}_${pin_short}.tsv" \
      assembly-region-reference "${repo_root}/${ref}" "${bam_resolved}" "${interval}" "${pad}"
  done <"${cases}"
fi

# B.5.1 region read payloads
cases="${repo_root}/parity/fixtures/hc-full-parity/b5-reads/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam interval padding _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    maybe_index_fa "${repo_root}/${ref}"
    bam_resolved="$(resolve_alignment_path "${repo_root}/${bam}" | tail -1)"
    maybe_index_bam "${bam_resolved}"
    pad="${padding:--}"
    write_dump "b5-reads/${case_id}_${pin_short}.tsv" \
      assembly-region-reads "${repo_root}/${ref}" "${bam_resolved}" "${interval}" "${pad}"
  done <"${cases}"
fi

# B.5.7 forceActive
cases="${repo_root}/parity/fixtures/hc-full-parity/b5-force/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam interval padding _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    maybe_index_fa "${repo_root}/${ref}"
    bam_resolved="$(resolve_alignment_path "${repo_root}/${bam}" | tail -1)"
    maybe_index_bam "${bam_resolved}"
    pad="${padding:--}"
    write_dump "b5-force/${case_id}_${pin_short}.tsv" \
      assembly-regions-force-active "${repo_root}/${ref}" "${bam_resolved}" "${interval}" "${pad}"
  done <"${cases}"
fi

# B.2 / B.3 / B.4 (optional padding column)
for phase in b2 b3 b4; do
  cases="${repo_root}/parity/fixtures/hc-full-parity/${phase}/cases.tsv"
  [[ -f "${cases}" ]] || continue
  subcmd=""
  case "${phase}" in
    b2) subcmd="assembly-regions" ;;
    b3) subcmd="apply-summary" ;;
    b4) subcmd="walker-traversal-summary" ;;
  esac
  while IFS=$'\t' read -r case_id ref bam interval padding _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    maybe_index_fa "${repo_root}/${ref}"
    bam_resolved="$(resolve_alignment_path "${repo_root}/${bam}" | tail -1)"
    maybe_index_bam "${bam_resolved}"
    pad="${padding:--}"
    write_dump "${phase}/${case_id}_${pin_short}.tsv" \
      "${subcmd}" "${repo_root}/${ref}" "${bam_resolved}" "${interval}" "${pad}"
  done <"${cases}"
done

# C.1 / C.2 / C.3 (no padding column in cases.tsv)
for phase in c1 c2 c3; do
  cases="${repo_root}/parity/fixtures/hc-full-parity/${phase}/cases.tsv"
  [[ -f "${cases}" ]] || continue
  subcmd=""
  case "${phase}" in
    c1) subcmd="raw-activity" ;;
    c2) subcmd="smoothed-activity" ;;
    c3) subcmd="active-locus" ;;
  esac
  while IFS=$'\t' read -r case_id ref bam interval _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    maybe_index_fa "${repo_root}/${ref}"
    bam_resolved="$(resolve_alignment_path "${repo_root}/${bam}" | tail -1)"
    maybe_index_bam "${bam_resolved}"
    write_dump "${phase}/${case_id}_${pin_short}.tsv" \
      "${subcmd}" "${repo_root}/${ref}" "${bam_resolved}" "${interval}" "-"
  done <"${cases}"
done

# C.4 genotype likelihoods (same intervals as C.1 subset)
cases="${repo_root}/parity/fixtures/hc-full-parity/c4-gl/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam interval _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    maybe_index_fa "${repo_root}/${ref}"
    bam_resolved="$(resolve_alignment_path "${repo_root}/${bam}" | tail -1)"
    maybe_index_bam "${bam_resolved}"
    write_dump "c4-gl/${case_id}_${pin_short}.tsv" \
      genotype-likelihoods "${repo_root}/${ref}" "${bam_resolved}" "${interval}" "-"
  done <"${cases}"
fi

# C.5 multisample joint raw-activity (same schema as C.1)
cases="${repo_root}/parity/fixtures/hc-full-parity/c5-multi/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam interval _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    maybe_index_fa "${repo_root}/${ref}"
    bam_resolved="$(resolve_alignment_path "${repo_root}/${bam}" | tail -1)"
    maybe_index_bam "${bam_resolved}"
    write_dump "c5-multi/${case_id}_${pin_short}.tsv" \
      raw-activity "${repo_root}/${ref}" "${bam_resolved}" "${interval}" "-"
  done <"${cases}"
fi

# Tier B C.5.3: force alleles + FeatureContext (raw-activity-force)
cases="${repo_root}/parity/fixtures/hc-full-parity/c5-force/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam interval vcf _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    maybe_index_fa "${repo_root}/${ref}"
    bam_resolved="$(resolve_alignment_path "${repo_root}/${bam}" | tail -1)"
    maybe_index_bam "${bam_resolved}"
    write_dump "c5-force/${case_id}_${pin_short}.tsv" \
      raw-activity-force "${repo_root}/${ref}" "${bam_resolved}" "${interval}" "${repo_root}/${vcf}" "-"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/c5-ploidy/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id sample_ploidy _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "c5-ploidy/${case_id}_${pin_short}.tsv" \
      ploidy-resolution "${sample_ploidy}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/g1/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id pairhmm_cases _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "g1/${case_id}_${pin_short}.tsv" \
      genotyping-aggregate "${repo_root}/${pairhmm_cases}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/g2-pl/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id fixture _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "g2-pl/${case_id}_${pin_short}.tsv" \
      genotype-format "${repo_root}/${fixture}"
  done <"${cases}"
fi

# H.1 reference confidence locus
cases="${repo_root}/parity/fixtures/hc-full-parity/h1/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam interval padding _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    maybe_index_fa "${repo_root}/${ref}"
    bam_resolved="$(resolve_alignment_path "${repo_root}/${bam}" | tail -1)"
    maybe_index_bam "${bam_resolved}"
    write_dump "h1/${case_id}_${pin_short}.tsv" \
      reference-confidence-locus \
      "${repo_root}/${ref}" "${bam_resolved}" "${interval}" "${padding}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/h1-inactive/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam interval padding _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    maybe_index_fa "${repo_root}/${ref}"
    bam_resolved="$(resolve_alignment_path "${repo_root}/${bam}" | tail -1)"
    maybe_index_bam "${bam_resolved}"
    write_dump "h1-inactive/${case_id}_${pin_short}.tsv" \
      inactive-reference-model \
      "${repo_root}/${ref}" "${bam_resolved}" "${interval}" "${padding}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/h2/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id contig length _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "h2/${case_id}_${pin_short}.tsv" gvcf-header "${contig}" "${length}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/h2-blocks/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id fixture _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "h2-blocks/${case_id}_${pin_short}.tsv" \
      gvcf-writer-blocks "${repo_root}/${fixture}"
  done <"${cases}"
fi

if [[ -f "${repo_root}/parity/fixtures/hc-full-parity/i1/expected/manifest.tsv" ]]; then
  write_dump "i1-manifest/manifest_${pin_short}.tsv" annotation-manifest
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/i1/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id alt_count samples _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "i1-core/${case_id}_${pin_short}.tsv" \
      annotate-core "${alt_count}" "${repo_root}/${samples}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/j2/cases_vcf.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id kind _ref bam_or_contig interval_or_pos ref_allele alt_allele gl_csv ad_csv _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    [[ "${kind}" == "call-region" ]] && continue
    write_dump "j2-vcf/${case_id}_${pin_short}.tsv" \
      variant-vcf-from-gl-ad \
      "${bam_or_contig}" "${interval_or_pos}" "${ref_allele}" "${alt_allele}" \
      "${gl_csv}" "${ad_csv}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/j2/cases_format.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id contig pos ref_allele alt_allele gl_csv ad_csv _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "j2-format/${case_id}_${pin_short}.tsv" \
      variant-format-from-gl-ad \
      "${contig}" "${pos}" "${ref_allele}" "${alt_allele}" "${gl_csv}" "${ad_csv}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/g2-region/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam interval padding target _active _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    maybe_index_fa "${repo_root}/${ref}"
    bam_resolved="$(resolve_alignment_path "${repo_root}/${bam}" | tail -1)"
    maybe_index_bam "${bam_resolved}"
    g2_args=(
      assembly-region-genotype
      "${repo_root}/${ref}" "${bam_resolved}" "${interval}"
    )
    [[ "${padding}" != "-" ]] && g2_args+=("${padding}")
    [[ -n "${target}" && "${target}" != "-" && "${target}" != "active" ]] && g2_args+=("${target}")
    write_dump "g2-region/${case_id}_${pin_short}.tsv" "${g2_args[@]}"
  done <"${cases}"
fi

# D.1
while IFS=$'\t' read -r case_id bam _expected; do
  [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
  bam_resolved="$(resolve_alignment_path "${repo_root}/${bam}" | tail -1)"
  write_dump "d1/${case_id}_${pin_short}.tsv" read-filters "${bam_resolved}"
done <"${repo_root}/parity/fixtures/hc-full-parity/d1/cases.tsv"

# D.2 positional downsamplers
cases="${repo_root}/parity/fixtures/hc-full-parity/d2/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam contig pos1 cap expected mode; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    bam_resolved="$(resolve_alignment_path "${repo_root}/${bam}" | tail -1)"
    maybe_index_bam "${bam_resolved}"
    if [[ -n "${mode:-}" ]]; then
      write_dump "d2/${case_id}_${pin_short}.tsv" \
        downsample-positional "${bam_resolved}" "${cap}" "${mode}"
    else
      write_dump "d2/${case_id}_${pin_short}.tsv" \
        downsample-positional "${bam_resolved}" "${cap}"
    fi
  done <"${cases}"
fi

# D.2.6 allele-biased downsampling
cases="${repo_root}/parity/fixtures/hc-full-parity/d2c/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id col2 col3 col4 col5 col6 col7; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    if [[ "${case_id}" == target_* ]]; then
      write_dump "d2c/${case_id}_${pin_short}.tsv" \
        allele-biased-target-counts "${col2}" "${col3}"
    else
      maybe_index_fa "${repo_root}/${col2}"
      bam_resolved="$(resolve_alignment_path "${repo_root}/${col3}" | tail -1)"
      maybe_index_bam "${bam_resolved}"
      write_dump "d2c/${case_id}_${pin_short}.tsv" \
        allele-biased-evidence "${repo_root}/${col2}" "${bam_resolved}" "${col4}" "${col5}" "${col6}"
    fi
  done <"${cases}"
fi

# D.2.7 contamination on isActive pileups
cases="${repo_root}/parity/fixtures/hc-full-parity/d2c-contam/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam interval contam _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    maybe_index_fa "${repo_root}/${ref}"
    bam_resolved="$(resolve_alignment_path "${repo_root}/${bam}" | tail -1)"
    maybe_index_bam "${bam_resolved}"
    write_dump "d2c-contam/${case_id}_${pin_short}.tsv" \
      raw-activity-contam "${repo_root}/${ref}" "${bam_resolved}" "${interval}" "${contam}"
  done <"${cases}"
fi

# D.3 HQ soft-clip mean (ReferenceConfidenceModel path)
cases="${repo_root}/parity/fixtures/hc-full-parity/d3/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam interval_cli _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    maybe_index_fa "${repo_root}/${ref}"
    bam_resolved="$(resolve_alignment_path "${repo_root}/${bam}" | tail -1)"
    maybe_index_bam "${bam_resolved}"
    write_dump "d3/${case_id}_${pin_short}.tsv" \
      soft-clip-mean "${repo_root}/${ref}" "${bam_resolved}" "${interval_cli}"
  done <"${cases}"
fi

# D.4 read shard pipeline (IUPAC pre / filter / post)
cases="${repo_root}/parity/fixtures/hc-full-parity/d4/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id bam _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    bam_resolved="$(resolve_alignment_path "${repo_root}/${bam}" | tail -1)"
    maybe_index_bam "${bam_resolved}"
    write_dump "d4/${case_id}_${pin_short}.tsv" read-shard-pipeline "${bam_resolved}"
  done <"${cases}"
fi

# PRE.1 soft-clip policy (revert vs hard-clip)
cases="${repo_root}/parity/fixtures/hc-full-parity/pre/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id bam _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    bam_resolved="$(resolve_alignment_path "${repo_root}/${bam}" | tail -1)"
    maybe_index_bam "${bam_resolved}"
    write_dump "pre/${case_id}_${pin_short}.tsv" read-pre-softclip "${bam_resolved}" 0 0
  done <"${cases}"
fi

# PRE.2 read length filter (unclippedReadLength >= 10)
cases="${repo_root}/parity/fixtures/hc-full-parity/pre-len/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id bam _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    bam_resolved="$(resolve_alignment_path "${repo_root}/${bam}" | tail -1)"
    maybe_index_bam "${bam_resolved}"
    write_dump "pre-len/${case_id}_${pin_short}.tsv" read-pre-len "${bam_resolved}"
  done <"${cases}"
fi

# PRE.3 assembly-path MQ filter (mapq >= 20)
cases="${repo_root}/parity/fixtures/hc-full-parity/pre-mq/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id bam _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    bam_resolved="$(resolve_alignment_path "${repo_root}/${bam}" | tail -1)"
    maybe_index_bam "${bam_resolved}"
    write_dump "pre-mq/${case_id}_${pin_short}.tsv" read-pre-mq "${bam_resolved}" 20
  done <"${cases}"
fi

# PRE.4 overlapping paired-fragment qual correction
cases="${repo_root}/parity/fixtures/hc-full-parity/pre-overlap/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id bam _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    bam_resolved="$(resolve_alignment_path "${repo_root}/${bam}" | tail -1)"
    maybe_index_bam "${bam_resolved}"
    write_dump "pre-overlap/${case_id}_${pin_short}.tsv" read-pre-overlap "${bam_resolved}"
  done <"${cases}"
fi

# E.1 assembly-graph edges
cases="${repo_root}/parity/fixtures/hc-full-parity/e1/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id reads kmer minq _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "e1/${case_id}_${pin_short}.tsv" \
      assembly-graph "${repo_root}/${reads}" "${kmer}" "${minq}"
  done <"${cases}"
fi

# E.2 multi-kmer edges
cases="${repo_root}/parity/fixtures/hc-full-parity/e2/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id reads kmer_csv minq _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "e2/${case_id}_${pin_short}.tsv" \
      assembly-graph-multi "${repo_root}/${reads}" "${kmer_csv}" "${minq}"
  done <"${cases}"
fi

# E.3 pruned graph summary
cases="${repo_root}/parity/fixtures/hc-full-parity/e3/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id reads kmer minq min_prune adaptive _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "e3/${case_id}_${pin_short}.tsv" \
      assembly-graph-summary "${repo_root}/${reads}" "${kmer}" "${minq}" "${min_prune}" "${adaptive}"
  done <"${cases}"
fi

# E.4 dangling recovery summary
cases="${repo_root}/parity/fixtures/hc-full-parity/e4/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref reads kmer minq min_prune min_dangling recover_heads _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "e4/${case_id}_${pin_short}.tsv" \
      assembly-graph-dangling-summary \
      "${repo_root}/${ref}" "${repo_root}/${reads}" \
      "${kmer}" "${minq}" "${min_prune}" "${min_dangling}" "${recover_heads}"
  done <"${cases}"
fi

# E.5 non-unique kmer / cycle policy summary
cases="${repo_root}/parity/fixtures/hc-full-parity/e5/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref reads kmer minq _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "e5/${case_id}_${pin_short}.tsv" \
      assembly-graph-non-unique-summary \
      "${ref}" "${repo_root}/${reads}" "${kmer}" "${minq}"
  done <"${cases}"
fi

# E.6 haplotype Smith–Waterman CIGARs
cases="${repo_root}/parity/fixtures/hc-full-parity/e6/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref haps _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "e6/${case_id}_${pin_short}.tsv" \
      assembly-haplotype-cigars \
      "${repo_root}/${ref}" "${repo_root}/${haps}"
  done <"${cases}"
fi

# E.7 assembler haplotype set (RT graph + GraphBasedKBest)
cases="${repo_root}/parity/fixtures/hc-full-parity/e7/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref reads kmer minq min_prune min_dangling recover_heads _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "e7/${case_id}_${pin_short}.tsv" \
      assembly-haplotypes \
      "${repo_root}/${ref}" "${repo_root}/${reads}" \
      "${kmer}" "${minq}" "${min_prune}" "${min_dangling}" "${recover_heads}"
  done <"${cases}"
fi

# E.7.1 KBest paths
cases="${repo_root}/parity/fixtures/hc-full-parity/e7-kbest/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref reads kmer minq min_prune min_dangling recover_heads max_haps _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "e7-kbest/${case_id}_${pin_short}.tsv" \
      assembly-kbest-paths \
      "${repo_root}/${ref}" "${repo_root}/${reads}" \
      "${kmer}" "${minq}" "${min_prune}" "${min_dangling}" "${recover_heads}" "${max_haps}"
  done <"${cases}"
fi

# E.7.2 haplotype cap
cases="${repo_root}/parity/fixtures/hc-full-parity/e7-cap/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref reads kmer minq min_prune min_dangling recover_heads max_haps _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "e7-cap/${case_id}_${pin_short}.tsv" \
      assembly-haplotypes-cap \
      "${repo_root}/${ref}" "${repo_root}/${reads}" \
      "${kmer}" "${minq}" "${min_prune}" "${min_dangling}" "${recover_heads}" "${max_haps}"
  done <"${cases}"
fi

# E.7.4 production ref tagging
cases="${repo_root}/parity/fixtures/hc-full-parity/e7-artificial/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref reads kmer minq min_prune min_dangling recover_heads _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "e7-artificial/${case_id}_${pin_short}.tsv" \
      assembly-haplotypes-production \
      "${repo_root}/${ref}" "${repo_root}/${reads}" \
      "${kmer}" "${minq}" "${min_prune}" "${min_dangling}" "${recover_heads}"
  done <"${cases}"
fi

# E2E assembly region haplotypes
cases="${repo_root}/parity/fixtures/hc-full-parity/e2e/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam interval padding target _l2 _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    maybe_index_fa "${repo_root}/${ref}"
    bam_resolved="$(resolve_alignment_path "${repo_root}/${bam}" | tail -1)"
    maybe_index_bam "${bam_resolved}"
    e2e_args=(
      assembly-region-haplotypes
      "${repo_root}/${ref}" "${bam_resolved}" "${interval}"
    )
    [[ "${padding}" != "-" ]] && e2e_args+=("${padding}")
    [[ -n "${target}" && "${target}" != "-" && "${target}" != "active" ]] && e2e_args+=("${target}")
    write_dump "e2e/${case_id}_${pin_short}.tsv" "${e2e_args[@]}"
  done <"${cases}"
fi

# E.0 full ReadThreadingAssembler multi-kmer assemble
cases="${repo_root}/parity/fixtures/hc-full-parity/e0-assemble/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref reads _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "e0-assemble/${case_id}_${pin_short}.tsv" \
      assembly-assemble \
      "${repo_root}/${ref}" "${repo_root}/${reads}"
  done <"${cases}"
fi

# E.1.1 pileup read error correction
cases="${repo_root}/parity/fixtures/hc-full-parity/e1-rec/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id reads log_odds _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "e1-rec/${case_id}_${pin_short}.tsv" \
      read-error-correction \
      "${repo_root}/${reads}" "${log_odds}"
  done <"${cases}"
fi

# E.8 SeqGraph post-processing summary
cases="${repo_root}/parity/fixtures/hc-full-parity/e8-seqgraph/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref reads kmer minq min_prune min_dangling recover_heads _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "e8-seqgraph/${case_id}_${pin_short}.tsv" \
      assembly-seqgraph-summary \
      "${repo_root}/${ref}" "${repo_root}/${reads}" \
      "${kmer}" "${minq}" "${min_prune}" "${min_dangling}" "${recover_heads}"
  done <"${cases}"
fi

# E.7.3 / E.7.5 JunctionTree KBest
for phase in e7-junction e7-edges; do
  cases="${repo_root}/parity/fixtures/hc-full-parity/${phase}/cases.tsv"
  if [[ ! -f "${cases}" ]]; then
    continue
  fi
  while IFS=$'\t' read -r case_id ref reads kmer minq recover_edges max_haps _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "${phase}/${case_id}_${pin_short}.tsv" \
      assembly-junction-haplotypes \
      "${repo_root}/${ref}" "${repo_root}/${reads}" \
      "${kmer}" "${minq}" "${recover_edges}" "${max_haps}"
  done <"${cases}"
done

# F.1 PairHMM likelihoods (scalar Java emitter aligned with Rust `pairhmm_log10_likelihood`).
cases="${repo_root}/parity/fixtures/hc-full-parity/f1/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id cases_path _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "f1/${case_id}_${pin_short}.tsv" \
      pairhmm-likelihoods "${repo_root}/${cases_path}"
  done <"${cases}"
fi

# F.2 GATK native Log10PairHMM (GKL when available).
cases="${repo_root}/parity/fixtures/hc-full-parity/f2/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id cases_path _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "f2/${case_id}_${pin_short}.tsv" \
      pairhmm-native-likelihoods "${repo_root}/${cases_path}"
  done <"${cases}"
fi

# E2E integration: region assembly → native PairHMM matrix.
cases="${repo_root}/parity/fixtures/hc-full-parity/e2e-int/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam interval padding target _l2 _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    maybe_index_fa "${repo_root}/${ref}"
    bam_resolved="$(resolve_alignment_path "${repo_root}/${bam}" | tail -1)"
    maybe_index_bam "${bam_resolved}"
    e2e_int_args=(
      assembly-region-pairhmm-likelihoods
      "${repo_root}/${ref}" "${bam_resolved}" "${interval}"
    )
    [[ "${padding}" != "-" ]] && e2e_int_args+=("${padding}")
    [[ -n "${target}" && "${target}" != "-" && "${target}" != "active" ]] && e2e_int_args+=("${target}")
    write_dump "e2e-int/${case_id}_${pin_short}.tsv" "${e2e_int_args[@]}"
  done <"${cases}"
fi

# F.3 BQ cap + haplotype filter gates.
cases="${repo_root}/parity/fixtures/hc-full-parity/f3/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r subphase case_id cases_path _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    case "${subphase}" in
      f3-bq)
        write_dump "f3/${case_id}_${pin_short}.tsv" \
          pairhmm-bq-cap "${repo_root}/${cases_path}"
        ;;
      f3-filter)
        write_dump "f3/${case_id}_${pin_short}.tsv" \
          pairhmm-haplotype-filter "${repo_root}/${cases_path}"
        ;;
    esac
  done <"${cases}"
fi

# --- Deferred gates (G/H/I/J/PRE L2) ---
g2af_tmp="${repo_root}/parity/reports/java-refresh-g2-af-tmp"
mkdir -p "${g2af_tmp}"
cases="${repo_root}/parity/fixtures/hc-full-parity/g2-af/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id gl_csv _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    input="${g2af_tmp}/${case_id}.input.tsv"
    printf '# case_id\tgl\n%s\t%s\n' "${case_id}" "${gl_csv}" >"${input}"
    write_dump "g2-af/${case_id}_${pin_short}.tsv" af-em "${input}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/g3/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ploidy max_gt _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "g3/${case_id}_${pin_short}.tsv" genotype-limits "${ploidy}" "${max_gt}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/g4/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id alleles phasing phase_set _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "g4/${case_id}_${pin_short}.tsv" genotype-phasing "${alleles}" "${phasing}" "${phase_set}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/g4-force/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id vcf contig pos filtered _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "g4-force/${case_id}_${pin_short}.tsv" \
      force-calling-genotype "${repo_root}/${vcf}" "${contig}" "${pos}" "${filtered}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/g-subset/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id sums is_ref max_alleles _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "g-subset/${case_id}_${pin_short}.tsv" \
      allele-subsetting "${sums}" "${is_ref}" "${max_alleles}"
  done <"${cases}"
fi

gsubpl_tmp="${repo_root}/parity/reports/java-refresh-g-subset-pl-tmp"
mkdir -p "${gsubpl_tmp}"
cases="${repo_root}/parity/fixtures/hc-full-parity/g-subset-pl/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id pl_csv ad_csv _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    input="${gsubpl_tmp}/${case_id}.input.tsv"
    printf '%s\t%s\t%s\n' "${case_id}" "${pl_csv}" "${ad_csv}" >"${input}"
    write_dump "g-subset-pl/${case_id}_${pin_short}.tsv" subset-alleles-pl "${input}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/g-subset-vc/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id _pl _ad _sac _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "g-subset-vc/${case_id}_${pin_short}.tsv" subset-alleles-vc \
      "${repo_root}/parity/fixtures/hc-full-parity/g-subset-vc/het_ac_sac_vc.tsv"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/g-subset-integration/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id sums is_ref max_alleles vc_fixture _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "g-subset-integration/${case_id}_${pin_short}.tsv" subset-alleles-integration \
      "${sums}" "${is_ref}" "${max_alleles}" "${repo_root}/${vc_fixture}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/h2-merge/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id merge_case _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "h2-merge/${case_id}_${pin_short}.tsv" \
      gvcf-merge-ref-confidence "${merge_case}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/g2-subset-live/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref bam interval padding target max_alleles assembly_profile _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    bam_resolved="$(resolve_alignment_path "${repo_root}/${bam}")"
    args=("${repo_root}/${ref}" "${bam_resolved}" "${interval}")
    if [[ -n "${padding}" && "${padding}" != "-" ]]; then
      args+=("${padding}")
    fi
    if [[ -n "${target}" && "${target}" != "-" ]]; then
      args+=("${target}")
    fi
    args+=("${assembly_profile:--}")
    args+=("${max_alleles}")
    write_dump "g2-subset-live/${case_id}_${pin_short}.tsv" \
      assembly-region-genotype-subset "${args[@]}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/h2-l5/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id fixture _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "h2-l5/${case_id}_${pin_short}.tsv" \
      gvcf-l5-merged "${repo_root}/${fixture}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/i1-standard/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref_fw ref_rv alt_fw alt_rv qual dp ref_bqs alt_bqs ref_pos alt_pos ref_mq alt_mq _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "i1-standard/${case_id}_${pin_short}.tsv" \
      standard-annotations "${ref_fw}" "${ref_rv}" "${alt_fw}" "${alt_rv}" "${qual}" "${dp}" \
      "${ref_bqs}" "${alt_bqs}" "${ref_pos}" "${alt_pos}" "${ref_mq}" "${alt_mq}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/i1-as/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id site_af site_qual _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "i1-as/${case_id}_${pin_short}.tsv" as-annotations "${site_af}" "${site_qual}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/i1-excess-het/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ref het hom _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "i1-excess-het/${case_id}_${pin_short}.tsv" excess-het "${ref}" "${het}" "${hom}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/i1-depth-hc/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id ad_csv _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "i1-depth-hc/${case_id}_${pin_short}.tsv" depth-per-sample-hc "${ad_csv}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/i1-plugins/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id plugin ref_fw ref_rv alt_fw alt_rv qual dp ref_bqs alt_bqs _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "i1-plugins/${case_id}_${pin_short}.tsv" \
      annotation-plugin "${plugin}" "${ref_fw}" "${ref_rv}" "${alt_fw}" "${alt_rv}" \
      "${qual}" "${dp}" "${ref_bqs}" "${alt_bqs}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/j-modes/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id mode has_variant locus_count _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "j-modes/${case_id}_${pin_short}.tsv" \
      emit-mode-decision "${mode}" "${has_variant}" "${locus_count}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/j-bamout/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id enabled write_count _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "j-bamout/${case_id}_${pin_short}.tsv" bamout-stub "${enabled}" "${write_count}"
  done <"${cases}"
fi

write_dump "j-dragen/default_off_${pin_short}.tsv" dragen-mode-branch

cases="${repo_root}/parity/fixtures/hc-full-parity/pre-dragstr/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id params_loaded _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "pre-dragstr/${case_id}_${pin_short}.tsv" dragstr-calibration "${params_loaded}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/e-debug/cases.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id failure_bam graph_dot _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    write_dump "e-debug/${case_id}_${pin_short}.tsv" assembly-debug-stub "${failure_bam}" "${graph_dot}"
  done <"${cases}"
fi

cases="${repo_root}/parity/fixtures/hc-full-parity/j2/cases_vcf.tsv"
if [[ -f "${cases}" ]]; then
  while IFS=$'\t' read -r case_id kind ref bam_or_contig interval_or_pos _ref_allele _alt_allele _gl _ad _expected; do
    [[ -z "${case_id}" || "${case_id}" == \#* ]] && continue
    if [[ "${kind}" == "call-region" ]]; then
      maybe_index_fa "${repo_root}/${ref}"
      bam_resolved="$(resolve_alignment_path "${repo_root}/${bam_or_contig}" | tail -1)"
      maybe_index_bam "${bam_resolved}"
      write_dump "j2-vcf/${case_id}_${pin_short}.tsv" \
        call-region-vcf "${repo_root}/${ref}" "${bam_resolved}" "${interval_or_pos}"
    fi
  done <"${cases}"
fi

echo "[hc-full-parity-java-refresh] wrote Java L2 dumps under ${out_root} (pin ${pin_short})"
