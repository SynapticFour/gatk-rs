#!/usr/bin/env bash
# Reclaim disk from Rust build artifacts and stale parity run logs.
#
# Safe to re-run. Does NOT touch:
#   - parity/reports/p12_diff/ (gate site lists)
#   - parity/reports/p12_realworld_na12878_20k.java.vcf (baseline)
#   - parity/reports/p12_l3_signoff_canonical.log
#   - parity/reports/p12_l4_signoff_canonical.log
#   - parity/reports/p12_l5_gvcf_canonical.log
#   - parity/reports/p12_l5_gvcf.json
#   - parity/reports/p12_l5_gvcf.md
#   - parity/reports/archive/
#   - parity/reports/hc-full-parity-l2/ (L2 signoff JSON)
#   - parity/realworld/ (BAM/FASTA assets)
#   - .gatk-src/ (local Java checkout; use --aggressive to remove)
# Usage: ./scripts/parity/clean_workspace.sh [--dry-run] [--aggressive]
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "${repo_root}"

dry_run=false
aggressive=false
for arg in "${1:-}" "${2:-}"; do
  case "${arg}" in
    --dry-run) dry_run=true; echo "DRY RUN — no deletes" ;;
    --aggressive) aggressive=true ;;
  esac
done

run_rm() {
  if $dry_run; then
    echo "would remove: $*"
  else
    rm -rf "$@"
  fi
}

echo "=== gatk-rs workspace cleanup ==="

# --- Rust build artifacts (~19 GB typical) ---
if [[ -d target ]]; then
  echo "cargo clean + remove target/"
  if $dry_run; then
    du -sh target 2>/dev/null || true
  else
    cargo clean 2>/dev/null || true
    rm -rf target
  fi
fi
run_rm target-agent target-foundation target-parity target-p10-release .docker-amd64-target

# External CARGO_TARGET_DIR used by parity scripts (~12 GB typical).
if [[ -d /tmp/gatk-parity-target ]]; then
  echo "remove /tmp/gatk-parity-target/"
  run_rm /tmp/gatk-parity-target
fi

if $aggressive && [[ -d .gatk-src ]]; then
  echo "aggressive: remove .gatk-src/ (re-clone via parity Java scripts)"
  run_rm .gatk-src
fi

# --- Stale parity/reports run logs (gitignored; safe to prune) ---
reports="${repo_root}/parity/reports"
if [[ -d "${reports}" ]]; then
  echo "prune stale logs under parity/reports/"

  # Canonical copies for docs (overwrite each cleanup)
  canonical_l3="${reports}/p12_l3_signoff_canonical.log"
  canonical_l4="${reports}/p12_l4_signoff_canonical.log"
  canonical_l5="${reports}/p12_l5_gvcf_canonical.log"
  canonical_trace="${reports}/p12_site_trace_latest.log"
  if ! $dry_run; then
    latest_l3="$(ls -t "${reports}"/p12_l3_signoff_20*.log 2>/dev/null | head -1 || true)"
    if [[ -n "${latest_l3}" && -f "${latest_l3}" ]]; then
      cp -f "${latest_l3}" "${canonical_l3}"
    fi
    latest_l4="$(ls -t "${reports}"/p12_l4_signoff_20*.log 2>/dev/null | head -1 || true)"
    if [[ -n "${latest_l4}" && -f "${latest_l4}" ]]; then
      cp -f "${latest_l4}" "${canonical_l4}"
    fi
    latest_l5="$(ls -t "${reports}"/p12_l5_gvcf_20*.log 2>/dev/null | head -1 || true)"
    if [[ -n "${latest_l5}" && -f "${latest_l5}" ]]; then
      cp -f "${latest_l5}" "${canonical_l5}"
    fi
    if [[ -f "${reports}/p12_java_site_trace_run2.log" ]]; then
      cp -f "${reports}/p12_java_site_trace_run2.log" "${canonical_trace}"
    fi
  fi

  # Drop timestamped sign-off copies once canonical exists.
  if ! $dry_run && [[ -f "${canonical_l3}" ]]; then
    for f in "${reports}"/p12_l3_signoff_20*.log; do
      [[ -e "${f}" ]] && run_rm "${f}"
    done
  fi
  if ! $dry_run && [[ -f "${canonical_l4}" ]]; then
    for f in "${reports}"/p12_l4_signoff_20*.log; do
      [[ -e "${f}" ]] && run_rm "${f}"
    done
  fi
  if ! $dry_run && [[ -f "${canonical_l5}" ]]; then
    for f in "${reports}"/p12_l5_gvcf_20*.log; do
      [[ -e "${f}" ]] && run_rm "${f}"
    done
  fi
  if ! $dry_run && [[ -f "${canonical_trace}" ]]; then
    for f in "${reports}"/p12_java_site_trace_run*.log; do
      [[ -e "${f}" ]] && run_rm "${f}"
    done
  fi

  # Intermediate ASM-8 / L3 / L4 iteration logs (canonical + key run logs kept).
  for f in "${reports}"/p12_asm8_production_run{5,6,7,8,9,10}.log \
           "${reports}"/p12_asm8_trace_short12.log \
           "${reports}"/p12_l3_phase3_run*.log \
           "${reports}"/p12_l3b_production_run{1,2,3,4}.log \
           "${reports}"/p12_l3a_regression_run9.log \
           "${reports}"/p12_java_site_trace_run1.log \
           "${reports}"/p12_six_sites_tmp.log \
           "${reports}"/p12_l3_signoff_20260609T112013Z.log \
           "${reports}"/p12_l3_signoff_20260609T122000Z.log \
           "${reports}"/p12_l3_signoff_20260609T122830Z.log \
           "${reports}"/p12_l3_v*.log \
           "${reports}"/p12_l3_l42*.log \
           "${reports}"/p12_l3_after_filter_fix.log \
           "${reports}"/p12_l3_gate_java_strict.log \
           "${reports}"/p12_l3_stabilization*.log \
           "${reports}"/p12_l3_site_trace_short*.log \
           "${reports}"/p12_l3_regression_hybrid.log \
           "${reports}"/p12_gate_66_v{2,3,4,5}.log \
           "${reports}"/p12_parity_gate_v*.log \
           "${reports}"/p12_parity_gate_final.log \
           "${reports}"/p12_92305634_*.log \
           "${reports}"/p12_site_92305634_*.log \
           "${reports}"/p12_92305716_*.log \
           "${reports}"/p12_site_92305716_*.log \
           "${reports}"/p12_two_sites*.log \
           "${reports}"/p12_emit_six_sites_after_java_strict.log \
           "${reports}"/p12_rust_only_triage.log \
           "${reports}"/p12_rust_only_triage_v2.log \
           "${reports}"/p12_format_parity_v*.log \
           "${reports}"/p12_format_parity_l42_v*.log \
           "${reports}"/p12_format_parity_l42_algorithmic.log \
           "${reports}"/p12_format_parity_l4fix.log \
           "${reports}"/p12_format_parity_l4.log \
           "${reports}"/p12_format_parity_hybrid.log \
           "${reports}"/p12_format_parity_fix*.log \
           "${reports}"/p12_format_parity_signoff.log \
           "${reports}"/p12_l3_gate_l4signoff.log \
           "${reports}"/p12_l3_gate_l4signoff_v2.log \
           "${reports}"/p12_l3_signoff_run.log \
           "${reports}"/p12_cluster_window_trace_extract.log \
           "${reports}"/p12_p13_realworld_full.log \
           "${reports}"/p12_gate_java_align*.log \
           "${reports}"/p12_l3_signoff_java_align.log \
           "${reports}"/p12_asm8_production_run11.log \
           "${reports}"/p12_l3a_regression_run10.log \
           "${reports}"/p12_l3b_production_run5.log \
           "${reports}"/p12_gate_66_latest.log \
           "${reports}"/p12_format_parity_signoff_v2.log \
           "${reports}"/p12_asm8_only_gate.log \
           "${reports}"/p12_format_parity_slice*.log \
           "${reports}"/p12_format_parity_16fix*.log \
           "${reports}"/p12_format_parity_tier3fix.log \
           "${reports}"/p12_format_parity_latest.log \
           "${reports}"/p12_format_parity_l42_run.log \
           "${reports}"/p12_l3_gate_l4signoff_v*.log; do
    [[ -e "${f}" ]] && run_rm "${f}"
  done

  # L2 tmp scratch + per-case JSON (regenerable; keep canonical log + summary).
  l2_dir="${reports}/hc-full-parity-l2"
  l2_tmp="${l2_dir}/tmp"
  if [[ -d "${l2_tmp}" ]]; then
    run_rm "${l2_tmp}"
  fi
  if [[ -d "${l2_dir}" ]]; then
    echo "prune L2 per-case JSON under hc-full-parity-l2/ (keep canonical log + l2_summary.json)"
    for f in "${l2_dir}"/*.json; do
      [[ -e "${f}" ]] || continue
      base="$(basename "${f}")"
      [[ "${base}" == "l2_summary.json" ]] && continue
      run_rm "${f}"
    done
  fi

  # Regenerable gVCF outputs from L5 battery.
  for f in "${reports}"/p12_l5_gvcf.java.g.vcf \
           "${reports}"/p12_l5_gvcf.java.g.vcf.idx \
           "${reports}"/p12_l5_gvcf.rust.g.vcf; do
    [[ -e "${f}" ]] && run_rm "${f}"
  done

  # Dev / live parity VCF scratch (keep realworld baselines for L6).
  for f in "${reports}"/live_*.vcf \
           "${reports}"/hc-realworld-*.vcf \
           "${reports}"/p11_*.vcf \
           "${reports}"/_p9_golden_gen.vcf; do
    [[ -e "${f}" ]] && run_rm "${f}"
  done

  # Parity script build scratch.
  run_rm "${repo_root}/parity/build"

  # Dev / CI scratch under parity/reports (regenerable).
  for d in realworld_pipeline_run p12_cluster_side_by_side asm_finalize_parity \
           hc-full-parity-java-refresh g-subset-pl-gen; do
    [[ -d "${reports}/${d}" ]] && run_rm "${reports}/${d}"
  done
  for pattern in countbases-* haplotypecaller-* printreads-* validate-help* \
                 hc-interval-list-* live_*.out live_*.live.* p3-* bam-alignment-*; do
    for f in "${reports}"/${pattern}; do
      [[ -e "${f}" ]] && run_rm "${f}"
    done
  done
  for f in "${reports}/hc-full-parity-l2_signoff.log" \
           "${reports}/p5_runtime_candidate_diff_details.json" \
           "${reports}/parity-smoke.json" \
           "${reports}/p6_pairhmm_live_drift_summary.json"; do
    [[ -e "${f}" ]] && run_rm "${f}"
  done

  # Stale analysis notes (canonical docs live under docs/ARCHITECTURE.md)
  for f in "${reports}"/L3_REGRESSION_*.md \
           "${reports}"/p12_l3_stabilization_site_trace_ANALYSIS.md \
           "${reports}"/p12_site_92305634_ANALYSIS.md; do
    [[ -e "${f}" ]] && run_rm "${f}"
  done

  # Parity script temp dirs
  for d in "${reports}"/hc-full-parity-*-tmp "${reports}"/java-refresh-*-tmp; do
    [[ -d "${d}" ]] && run_rm "${d}"
  done
fi

echo "=== cleanup done ==="
if ! $dry_run; then
  du -sh target 2>/dev/null || echo "target/ removed"
  du -sh .gatk-src 2>/dev/null || echo ".gatk-src/ removed (re-fetch via parity Java scripts)"
  du -sh parity/reports 2>/dev/null || true
  du -sh . 2>/dev/null || true
fi
