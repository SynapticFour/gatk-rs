#!/usr/bin/env bash
# Nightly GIAB Ashkenazi trio equivalence (Prompt-B spine):
#   HaplotypeCaller (-ERC GVCF) → CombineGVCFs → GenotypeGVCFs → VariantFiltration
# on chr20 + chr21 + hard GIAB-stratification slices.
#
# Disk-safe: never downloads full WGS BAMs. For each region and sample:
#   curl BAI → samtools view -L <bed> -X <remote.bam> <bai> → HC GVCF → rm BAM
#
# Env (selected):
#   NIGHTLY_OUT_ROOT, NIGHTLY_HARD_BUDGET_BP, NIGHTLY_SKIP_JAVA
#   NIGHTLY_STAND_CALL_CONF, HAPPY_BIN / HAPPY_DOCKER_IMAGE
#   GIAB_* BAM URL overrides, CARGO_TARGET_DIR
set -euo pipefail

repo_root="$(cd "$(dirname "$0")/../../.." && pwd)"
cd "${repo_root}"
# shellcheck source=../lib_pinned_gatk.sh
source "${repo_root}/scripts/parity/lib_pinned_gatk.sh"
# shellcheck source=../m4_disk_guard.sh
source "${repo_root}/scripts/parity/m4_disk_guard.sh"

stamp="$(date -u +%Y%m%dT%H%M%SZ)"
out_root="${NIGHTLY_OUT_ROOT:-${repo_root}/parity/giab/runs/nightly_trio_${stamp}}"
truth_root="${GIAB_TRUTH_ROOT:-${repo_root}/parity/giab/truth}"
strat_root="${GIAB_STRAT_ROOT:-${repo_root}/parity/giab/stratifications/GRCh37}"
ref="${GIAB_REFERENCE:-${repo_root}/parity/realworld/assets/hs37d5.simple.fa}"
regions_dir="${out_root}/region_defs"
hard_budget="${NIGHTLY_HARD_BUDGET_BP:-2000000}"
stand_call_conf="${NIGHTLY_STAND_CALL_CONF:-30}"
skip_java="${NIGHTLY_SKIP_JAVA:-0}"
happy_image="${HAPPY_DOCKER_IMAGE:-jmcdani20/hap.py:v0.3.12}"
happy_entry="${HAPPY_DOCKER_ENTRYPOINT:-/opt/hap.py/bin/hap.py}"
ftp_data="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data"

mkdir -p "${out_root}"
log="${out_root}/run.log"
exec > >(tee -a "${log}") 2>&1

echo "=== nightly trio equivalence ${stamp} ==="
echo "out_root=${out_root}"
df -h . || true

target_dir="${CARGO_TARGET_DIR:-${repo_root}/target}"
rust_bin="${NIGHTLY_RUST_BIN:-}"
if [[ -z "${rust_bin}" ]]; then
  for cand in "${target_dir}/release/gatk-rs" "${target_dir}/debug/gatk-rs"; do
    [[ -x "${cand}" ]] && rust_bin="${cand}" && break
  done
fi
if [[ -z "${rust_bin}" || ! -x "${rust_bin}" ]]; then
  echo "[nightly-trio] building gatk-rs release…"
  "${repo_root}/scripts/parity/build_gatk_rs_release.sh"
  rust_bin="${target_dir}/release/gatk-rs"
fi

if ! command -v samtools >/dev/null 2>&1; then
  echo "[nightly-trio] samtools required" >&2
  exit 2
fi
if ! command -v docker >/dev/null 2>&1 && [[ -z "${HAPPY_BIN:-}" ]]; then
  echo "[nightly-trio] docker or HAPPY_BIN required for hap.py" >&2
  exit 2
fi

# --- Truth / strat / reference ------------------------------------------------
if [[ "${GIAB_FETCH_TRUTH:-1}" == "1" ]]; then
  "${repo_root}/scripts/parity/giab/fetch_giab_truthsets.sh"
fi
if [[ "${GIAB_STAGE_REF:-1}" == "1" && ! -f "${ref}" ]]; then
  "${repo_root}/scripts/parity/realworld/03_stage_reference_and_truth.sh"
  ref="${repo_root}/parity/realworld/assets/hs37d5.simple.fa"
fi
[[ -f "${ref}" ]] || { echo "missing ref ${ref}" >&2; exit 2; }

truth_vcf="${truth_root}/HG002_GRCh37_1_22_v4.2.1_benchmark.vcf.gz"
truth_bed="${truth_root}/HG002_GRCh37_1_22_v4.2.1_benchmark.bed"
[[ -f "${truth_vcf}" && -f "${truth_bed}" ]] || {
  echo "missing HG002 truth under ${truth_root}" >&2
  exit 2
}

python3 "${repo_root}/scripts/parity/giab/build_nightly_regions.py" \
  --strat-root "${strat_root}" \
  --out-dir "${regions_dir}" \
  --hard-budget-bp "${hard_budget}"

# --- BAM URLs (AJ trio 300x novoalign / hs37d5) --------------------------------
bam_url_hg002="${GIAB_HG002_BAM_URL:-${ftp_data}/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.300x.bam}"
bam_url_hg003="${GIAB_HG003_BAM_URL:-${ftp_data}/AshkenazimTrio/HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG003.hs37d5.300x.bam}"
bam_url_hg004="${GIAB_HG004_BAM_URL:-${ftp_data}/AshkenazimTrio/HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG004.hs37d5.300x.bam}"

samples=(HG002 HG003 HG004)
bam_urls=("${bam_url_hg002}" "${bam_url_hg003}" "${bam_url_hg004}")

# --- hap.py wrapper -----------------------------------------------------------
bin_dir="${out_root}/bin"
mkdir -p "${bin_dir}"
happy_wrap="${bin_dir}/hap.py"
if [[ -n "${HAPPY_BIN:-}" ]]; then
  ln -sfn "${HAPPY_BIN}" "${happy_wrap}"
else
  cat >"${happy_wrap}" <<EOF
#!/usr/bin/env bash
set -euo pipefail
exec docker run --rm --platform linux/amd64 \\
  -v "${repo_root}:${repo_root}" \\
  -w "${repo_root}" \\
  "${happy_image}" \\
  "${happy_entry}" "\$@"
EOF
  chmod +x "${happy_wrap}"
  if [[ "${NIGHTLY_SKIP_DOCKER_PULL:-0}" != "1" ]]; then
    docker pull --platform linux/amd64 "${happy_image}" || true
    docker pull --platform "${GATK_DOCKER_PLATFORM}" "${GATK_DOCKER_IMAGE}" || true
  fi
fi

slice_bam() {
  local sample="$1" bam_url="$2" bed="$3" out_bam="$4"
  local bai_local="${out_bam}.remote.bai"
  mkdir -p "$(dirname "${out_bam}")"
  if [[ ! -f "${bai_local}" ]]; then
    echo "[nightly-trio] BAI ${sample}…"
    curl -fL --retry 3 --retry-delay 2 -o "${bai_local}.partial" "${bam_url}.bai"
    mv -f "${bai_local}.partial" "${bai_local}"
  fi
  echo "[nightly-trio] samtools view -L $(basename "${bed}") ${sample}…"
  # -L: BED regions only; -X: remote BAM + local BAI (no full WGS download)
  samtools view -b -L "${bed}" -X "${bam_url}" "${bai_local}" -o "${out_bam}.partial"
  mv -f "${out_bam}.partial" "${out_bam}"
  samtools index "${out_bam}"
}

extract_sample_vcf() {
  local in_vcf="$1" sample="$2" out_vcf="$3"
  # Prefer bcftools; fall back to a tiny awk that keeps the named sample column.
  if command -v bcftools >/dev/null 2>&1; then
    bcftools view -s "${sample}" -O v -o "${out_vcf}" "${in_vcf}"
    return 0
  fi
  python3 - "${in_vcf}" "${sample}" "${out_vcf}" <<'PY'
import sys
inp, sample, outp = sys.argv[1:4]
with open(inp, encoding="utf-8", errors="replace") as fh, open(outp, "w", encoding="utf-8") as out:
    sample_idx = None
    for line in fh:
        if line.startswith("##"):
            out.write(line)
            continue
        if line.startswith("#CHROM"):
            cols = line.rstrip("\n").split("\t")
            try:
                sample_idx = cols.index(sample)
            except ValueError:
                # GIAB BAMs often use SM=HG002; joint VCF may use that or a path stem.
                for i, c in enumerate(cols[9:], start=9):
                    if sample in c:
                        sample_idx = i
                        break
            if sample_idx is None:
                raise SystemExit(f"sample {sample} not in VCF header: {cols[9:]}")
            out.write("\t".join(cols[:9] + [sample]) + "\n")
            continue
        cols = line.rstrip("\n").split("\t")
        out.write("\t".join(cols[:9] + [cols[sample_idx]]) + "\n")
PY
}

run_region() {
  local name="$1" bed="$2" intervals_file="$3" kind="$4" span_bp="$5"
  local rdir="${out_root}/regions/${name}"
  mkdir -p "${rdir}"
  cat >"${rdir}/region.json" <<EOF
{"name": "${name}", "kind": "${kind}", "span_bp": ${span_bp}, "bed": "${bed}", "status": "running"}
EOF

  echo "======== region ${name} (${kind}, ${span_bp} bp) ========"
  m4_require_free_gb "${NIGHTLY_MIN_FREE_GB:-6}" || {
    echo "[nightly-trio] low disk before ${name}; continuing best-effort"
  }

  # Interval arg for callers: prefer intervals file (1-based); GATK/rust accept path.
  local interval_arg="${intervals_file}"

  local rust_gvcfs=()
  local java_gvcfs=()
  local i
  for i in 0 1 2; do
    local sample="${samples[$i]}"
    local bam_url="${bam_urls[$i]}"
    local bam="${rdir}/bam/${sample}.bam"
    slice_bam "${sample}" "${bam_url}" "${bed}" "${bam}"

    local rg="${rdir}/rust.${sample}.g.vcf"
    echo "[nightly-trio] Rust HC GVCF ${sample} @ ${name}"
    "${rust_bin}" haplotypecaller \
      -R "${ref}" -I "${bam}" -O "${rg}" -L "${interval_arg}" \
      --emit-ref-confidence GVCF
    rust_gvcfs+=("${rg}")

    if [[ "${skip_java}" != "1" ]]; then
      local jg="${rdir}/java.${sample}.g.vcf"
      echo "[nightly-trio] Java HC GVCF ${sample} @ ${name}"
      "${repo_root}/scripts/parity/run_java_gatk.sh" \
        "${rdir}/java.hc.${sample}.stdout" \
        HaplotypeCaller -R "${ref}" -I "${bam}" -O "${jg}" -L "${interval_arg}" \
        -ERC GVCF --verbosity ERROR
      "${repo_root}/scripts/parity/run_java_gatk.sh" \
        "${rdir}/java.idx.${sample}.stdout" \
        IndexFeatureFile -I "${jg}" || true
      java_gvcfs+=("${jg}")
    fi

    # Free BAM slice immediately (largest disk consumer)
    rm -f "${bam}" "${bam}.bai"
  done

  local rust_combined="${rdir}/rust.combined.g.vcf"
  local rust_gt="${rdir}/rust.genotyped.vcf"
  local rust_filt="${rdir}/rust.filtered.vcf"
  "${rust_bin}" combine-gvcfs -R "${ref}" -L "${interval_arg}" \
    -V "${rust_gvcfs[0]}" -V "${rust_gvcfs[1]}" -V "${rust_gvcfs[2]}" \
    -O "${rust_combined}"
  "${rust_bin}" genotype-gvcfs -R "${ref}" -V "${rust_combined}" -O "${rust_gt}" \
    --stand-call-conf "${stand_call_conf}"
  "${rust_bin}" variant-filtration -V "${rust_gt}" -O "${rust_filt}" --preset snp

  if [[ "${skip_java}" != "1" ]]; then
    local java_combined="${rdir}/java.combined.g.vcf"
    local java_gt="${rdir}/java.genotyped.vcf"
    local java_filt="${rdir}/java.filtered.vcf"
    "${repo_root}/scripts/parity/run_java_gatk.sh" \
      "${rdir}/java.combine.stdout" \
      CombineGVCFs -R "${ref}" -L "${interval_arg}" \
      -V "${java_gvcfs[0]}" -V "${java_gvcfs[1]}" -V "${java_gvcfs[2]}" \
      -O "${java_combined}"
    "${repo_root}/scripts/parity/run_java_gatk.sh" \
      "${rdir}/java.gg.stdout" \
      GenotypeGVCFs -R "${ref}" -V "${java_combined}" -O "${java_gt}" \
      --standard-min-confidence-threshold-for-calling "${stand_call_conf}"
    "${repo_root}/scripts/parity/run_java_gatk.sh" \
      "${rdir}/java.vf.stdout" \
      VariantFiltration -V "${java_gt}" -O "${java_filt}" \
      --filter-expression "QD < 2.0" --filter-name "QD2" \
      --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
      --filter-expression "SOR > 3.0" --filter-name "SOR3" \
      --filter-expression "FS > 60.0" --filter-name "FS60" \
      --filter-expression "MQ < 40.0" --filter-name "MQ40" \
      --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
      --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8"
  fi

  # Score HG002 column vs GIAB HG002 truth (intersection with region BED ∩ confident)
  local conf_region="${rdir}/eval_confident.bed"
  # Intersect truth confident BED with calling BED (python; bedtools optional)
  python3 - "${truth_bed}" "${bed}" "${conf_region}" <<'PY'
import sys
from pathlib import Path

def load(p):
    rows = []
    for line in Path(p).read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.strip() or line.startswith("#") or line.startswith("track"):
            continue
        c = line.split("\t")
        if len(c) < 3:
            continue
        chrom = c[0][3:] if c[0].startswith("chr") else c[0]
        rows.append((chrom, int(c[1]), int(c[2])))
    return rows

def intersect(a, b):
    b_by = {}
    for chrom, s, e in b:
        b_by.setdefault(chrom, []).append((s, e))
    for chrom in b_by:
        b_by[chrom].sort()
    out = []
    for chrom, s, e in a:
        for bs, be in b_by.get(chrom, []):
            if be <= s:
                continue
            if bs >= e:
                break
            out.append((chrom, max(s, bs), min(e, be)))
    return out

truth, call, dest = sys.argv[1:4]
iv = intersect(load(truth), load(call))
Path(dest).write_text("".join(f"{c}\t{s}\t{e}\n" for c, s, e in iv), encoding="utf-8")
if not iv:
    # Fall back to calling BED so hap.py still runs
    Path(dest).write_text(Path(call).read_text(encoding="utf-8"), encoding="utf-8")
PY

  local rust_hg002="${rdir}/rust.HG002.vcf"
  extract_sample_vcf "${rust_filt}" "HG002" "${rust_hg002}" || \
    extract_sample_vcf "${rust_filt}" "sample" "${rust_hg002}" || true
  if [[ ! -s "${rust_hg002}" ]]; then
    # Last resort: use multi-sample filtered VCF
    cp -f "${rust_filt}" "${rust_hg002}"
  fi

  mkdir -p "${rdir}/happy_rust"
  set +e
  "${happy_wrap}" "${truth_vcf}" "${rust_hg002}" \
    -r "${ref}" -f "${conf_region}" \
    -o "${rdir}/happy_rust/prefix" \
    --threads "${NIGHTLY_HAPPY_THREADS:-2}"
  local happy_rc=$?
  set -e

  if [[ "${skip_java}" != "1" && -f "${rdir}/java.filtered.vcf" ]]; then
    local java_hg002="${rdir}/java.HG002.vcf"
    extract_sample_vcf "${rdir}/java.filtered.vcf" "HG002" "${java_hg002}" || \
      cp -f "${rdir}/java.filtered.vcf" "${java_hg002}"
    mkdir -p "${rdir}/happy_java"
    set +e
    "${happy_wrap}" "${truth_vcf}" "${java_hg002}" \
      -r "${ref}" -f "${conf_region}" \
      -o "${rdir}/happy_java/prefix" \
      --threads "${NIGHTLY_HAPPY_THREADS:-2}"
    set -e
  fi

  # Drop bulky intermediates; keep filtered VCFs + hap.py outputs
  rm -f "${rdir}/rust."*.g.vcf "${rdir}/java."*.g.vcf \
        "${rdir}/rust.combined.g.vcf" "${rdir}/java.combined.g.vcf" || true

  local status="ok"
  if [[ "${happy_rc}" -ne 0 || ! -f "${rdir}/happy_rust/prefix.summary.csv" ]]; then
    status="happy_failed"
  fi
  python3 - "${rdir}/region.json" "${status}" <<'PY'
import json, pathlib, sys
p, status = pathlib.Path(sys.argv[1]), sys.argv[2]
data = json.loads(p.read_text(encoding="utf-8"))
data["status"] = status
p.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
PY
  echo "[nightly-trio] region ${name} status=${status}"
}

# --- Drive all regions (best-effort; soft failures) ---------------------------
overall_rc=0
while IFS= read -r region_json; do
  [[ -n "${region_json}" ]] || continue
  name="$(printf '%s' "${region_json}" | python3 -c 'import json,sys; print(json.load(sys.stdin)["name"])')"
  bed="$(printf '%s' "${region_json}" | python3 -c 'import json,sys; print(json.load(sys.stdin)["bed"])')"
  ivf="$(printf '%s' "${region_json}" | python3 -c 'import json,sys; print(json.load(sys.stdin)["intervals_file"])')"
  kind="$(printf '%s' "${region_json}" | python3 -c 'import json,sys; print(json.load(sys.stdin)["kind"])')"
  span="$(printf '%s' "${region_json}" | python3 -c 'import json,sys; print(json.load(sys.stdin)["span_bp"])')"
  set +e
  run_region "${name}" "${bed}" "${ivf}" "${kind}" "${span}"
  rc=$?
  set -e
  if [[ "${rc}" -ne 0 ]]; then
    echo "[nightly-trio] WARNING: region ${name} failed rc=${rc}"
    overall_rc=1
    mkdir -p "${out_root}/regions/${name}"
    echo "{\"name\": \"${name}\", \"kind\": \"${kind}\", \"span_bp\": ${span}, \"status\": \"failed\"}" \
      >"${out_root}/regions/${name}/region.json"
  fi
done < <(python3 -c 'import json,pathlib,sys; m=json.loads(pathlib.Path(sys.argv[1]).read_text(encoding="utf-8"));
[print(json.dumps(r)) for r in m["regions"]]' "${regions_dir}/manifest.json")

# Publish summary hooks are invoked by the workflow; still write a local pointer.
echo "${out_root}" >"${repo_root}/parity/giab/runs/nightly_trio_latest.txt"
echo "[nightly-trio] done out_root=${out_root} overall_rc=${overall_rc}"
# Soft gate: always exit 0 from the orchestrator when regions were attempted;
# regression issues are handled after publish. Non-zero only if zero regions ran.
if [[ ! -d "${out_root}/regions" ]]; then
  exit 2
fi
exit 0
