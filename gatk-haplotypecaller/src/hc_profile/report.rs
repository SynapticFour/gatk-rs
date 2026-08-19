//! JSON + Markdown report writers for HC production profiling.

use std::fmt::Write as FmtWrite;
use std::fs;
use std::path::Path;
use std::time::Duration;

use super::genotype::GenotypeAgg;
use super::pairhmm::PairHmmAgg;
use super::stages::{Stage, StageAgg};

pub fn write_reports(
    json_path: &Path,
    stages: &StageAgg,
    pairhmm: &PairHmmAgg,
    genotype: &GenotypeAgg,
    run_wall: Duration,
    run_cpu: Option<Duration>,
    rayon_threads: usize,
    regions: u64,
) -> std::io::Result<()> {
    if let Some(parent) = json_path.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)?;
        }
    }
    let json = build_json(
        stages,
        pairhmm,
        genotype,
        run_wall,
        run_cpu,
        rayon_threads,
        regions,
    );
    fs::write(json_path, &json)?;
    let md_path = json_path.with_extension("md");
    let md = build_markdown(
        stages,
        pairhmm,
        genotype,
        run_wall,
        run_cpu,
        rayon_threads,
        regions,
    );
    fs::write(md_path, md)?;
    Ok(())
}

fn ns_to_s(ns: u64) -> f64 {
    ns as f64 / 1e9
}

fn build_json(
    stages: &StageAgg,
    pairhmm: &PairHmmAgg,
    genotype: &GenotypeAgg,
    run_wall: Duration,
    run_cpu: Option<Duration>,
    rayon_threads: usize,
    regions: u64,
) -> String {
    let run_wall_s = run_wall.as_secs_f64();
    let run_cpu_s = run_cpu.map(|d| d.as_secs_f64());
    let parallel_eff = match run_cpu_s {
        Some(cpu) if run_wall_s > 0.0 => Some(cpu / (run_wall_s * rayon_threads as f64)),
        _ => None,
    };
    let occ = pairhmm.simd_occupancy();

    let mut out = String::new();
    out.push_str("{\n");
    out.push_str("  \"schema\": \"gatk-rs.hc_profile.v1\",\n");
    out.push_str("  \"honesty\": \"Observe-only production profile. Does not change genotype/emit. Stage walls may overlap TRACE dual-write; prefer PairHMM/genotype nested counters for leaf detail.\",\n");
    out.push_str(&format!("  \"regions\": {regions},\n"));
    out.push_str(&format!("  \"rayon_threads\": {rayon_threads},\n"));
    out.push_str(&format!("  \"run_wall_s\": {run_wall_s:.6},\n"));
    match run_cpu_s {
        Some(c) => out.push_str(&format!("  \"run_cpu_s\": {c:.6},\n")),
        None => out.push_str("  \"run_cpu_s\": null,\n"),
    }
    match parallel_eff {
        Some(e) => out.push_str(&format!("  \"parallel_efficiency\": {e:.4},\n")),
        None => out.push_str("  \"parallel_efficiency\": null,\n"),
    }
    out.push_str("  \"stages\": {\n");
    for (i, stage) in Stage::all().iter().enumerate() {
        let s = stages.get(*stage);
        let comma = if i + 1 == Stage::all().len() { "" } else { "," };
        out.push_str(&format!(
            "    \"{}\": {{\n      \"calls\": {},\n      \"wall_s\": {:.6},\n      \"cpu_s\": {:.6},\n      \"cpu_samples\": {},\n      \"avg_wall_s\": {:.9},\n      \"alloc_bytes\": {},\n      \"alloc_events\": {}\n    }}{comma}\n",
            stage.as_str(),
            s.calls,
            ns_to_s(s.wall_ns),
            ns_to_s(s.cpu_ns),
            s.cpu_samples,
            s.avg_wall_ns() / 1e9,
            s.alloc_bytes,
            s.alloc_events
        ));
    }
    out.push_str("  },\n");
    out.push_str("  \"pairhmm\": {\n");
    out.push_str(&format!("    \"calls\": {},\n", pairhmm.calls));
    out.push_str(&format!(
        "    \"reads_scored\": {},\n",
        pairhmm.reads_scored
    ));
    out.push_str(&format!(
        "    \"haplotypes_scored\": {},\n",
        pairhmm.haplotypes_scored
    ));
    out.push_str(&format!(
        "    \"read_haplotype_pairs\": {},\n",
        pairhmm.read_hap_pairs
    ));
    out.push_str(&format!(
        "    \"mean_haplotypes_per_read\": {:.4},\n",
        pairhmm.mean_haps_per_read()
    ));
    out.push_str(&format!(
        "    \"mean_read_len\": {:.2},\n",
        if pairhmm.reads_scored == 0 {
            0.0
        } else {
            pairhmm.read_len_sum as f64 / pairhmm.reads_scored as f64
        }
    ));
    out.push_str(&format!(
        "    \"mean_hap_len\": {:.2},\n",
        if pairhmm.haplotypes_scored == 0 {
            0.0
        } else {
            pairhmm.hap_len_sum as f64 / pairhmm.haplotypes_scored as f64
        }
    ));
    out.push_str(&format!("    \"simd_pack_units\": {},\n", occ.pack_units));
    out.push_str(&format!("    \"simd_pack_lanes\": {},\n", occ.pack_lanes));
    out.push_str(&format!(
        "    \"simd_pack_haps_est\": {},\n",
        occ.pack_hap_est
    ));
    out.push_str(&format!(
        "    \"simd_pack_occupancy_pct\": {:.2},\n",
        occ.pack_occupancy_pct
    ));
    out.push_str(&format!(
        "    \"prefix_reuse_haps\": {},\n",
        pairhmm.prefix_reuse_haps
    ));
    out.push_str(&format!(
        "    \"leftover_haps\": {},\n",
        pairhmm.leftover_haps
    ));
    out.push_str(&format!(
        "    \"prefix_reuse_pct\": {:.2},\n",
        occ.prefix_reuse_pct
    ));
    out.push_str(&format!("    \"leftover_pct\": {:.2},\n", occ.leftover_pct));
    out.push_str(&format!(
        "    \"dp_cells_evaluated\": {},\n",
        pairhmm.dp_cells_evaluated
    ));
    out.push_str(&format!(
        "    \"dp_cells_avoided_prefix\": {},\n",
        pairhmm.dp_cells_avoided_prefix
    ));
    out.push_str(&format!(
        "    \"wall_s\": {:.6},\n",
        ns_to_s(pairhmm.wall_ns)
    ));
    out.push_str("    \"read_len_hist_bucket25\": {\n");
    write_hist(&mut out, &pairhmm.read_len_hist);
    out.push_str("    },\n");
    out.push_str("    \"hap_len_hist_bucket25\": {\n");
    write_hist(&mut out, &pairhmm.hap_len_hist);
    out.push_str("    }\n");
    out.push_str("  },\n");
    out.push_str("  \"genotype\": {\n");
    out.push_str(&format!("    \"sites\": {},\n", genotype.sites));
    out.push_str(&format!(
        "    \"candidate_alleles_sum\": {},\n",
        genotype.candidate_alleles_sum
    ));
    out.push_str(&format!(
        "    \"genotype_states_sum\": {},\n",
        genotype.genotype_states_sum
    ));
    out.push_str(&format!(
        "    \"pl_vector_len_sum\": {},\n",
        genotype.pl_vector_len_sum
    ));
    out.push_str(&format!("    \"samples_sum\": {},\n", genotype.samples_sum));
    out.push_str(&format!(
        "    \"mean_alleles_per_site\": {:.4},\n",
        if genotype.sites == 0 {
            0.0
        } else {
            genotype.candidate_alleles_sum as f64 / genotype.sites as f64
        }
    ));
    out.push_str(&format!(
        "    \"mean_states_per_site\": {:.4},\n",
        if genotype.sites == 0 {
            0.0
        } else {
            genotype.genotype_states_sum as f64 / genotype.sites as f64
        }
    ));
    out.push_str(&format!(
        "    \"time_per_site_s\": {:.9},\n",
        genotype.time_per_site_ns() / 1e9
    ));
    out.push_str(&format!(
        "    \"time_per_genotype_state_s\": {:.9},\n",
        genotype.time_per_state_ns() / 1e9
    ));
    out.push_str(&format!(
        "    \"wall_s\": {:.6},\n",
        ns_to_s(genotype.wall_ns)
    ));
    out.push_str(&format!(
        "    \"ad_wall_s\": {:.6},\n",
        ns_to_s(genotype.ad_wall_ns)
    ));
    out.push_str(&format!(
        "    \"allele_map_wall_s\": {:.6},\n",
        ns_to_s(genotype.allele_map_wall_ns)
    ));
    out.push_str(&format!(
        "    \"marginalize_wall_s\": {:.6},\n",
        ns_to_s(genotype.marginalize_wall_ns)
    ));
    out.push_str(&format!(
        "    \"genotype_enum_wall_s\": {:.6},\n",
        ns_to_s(genotype.genotype_enum_wall_ns)
    ));
    out.push_str(&format!(
        "    \"event_rebuild_wall_s\": {:.6}\n",
        ns_to_s(genotype.event_rebuild_wall_ns)
    ));
    out.push_str("  }\n");
    out.push_str("}\n");
    out
}

fn write_hist(out: &mut String, hist: &std::collections::BTreeMap<u32, u64>) {
    let items: Vec<_> = hist.iter().collect();
    for (i, (k, v)) in items.iter().enumerate() {
        let comma = if i + 1 == items.len() { "" } else { "," };
        let _ = writeln!(out, "      \"{k}\": {v}{comma}");
    }
}

fn build_markdown(
    stages: &StageAgg,
    pairhmm: &PairHmmAgg,
    genotype: &GenotypeAgg,
    run_wall: Duration,
    run_cpu: Option<Duration>,
    rayon_threads: usize,
    regions: u64,
) -> String {
    let mut md = String::new();
    md.push_str("# HC production profile\n\n");
    md.push_str("Observe-only. Does **not** change genotype/emit.\n\n");
    md.push_str(&format!("- Regions: **{regions}**\n"));
    md.push_str(&format!("- Rayon threads: **{rayon_threads}**\n"));
    md.push_str(&format!("- Run wall: **{:.3}s**\n", run_wall.as_secs_f64()));
    if let Some(cpu) = run_cpu {
        let eff = cpu.as_secs_f64() / (run_wall.as_secs_f64() * rayon_threads as f64).max(1e-9);
        md.push_str(&format!(
            "- Run CPU (user+sys): **{:.3}s** (parallel efficiency ≈ **{:.2}**)\n",
            cpu.as_secs_f64(),
            eff
        ));
    }
    md.push_str("\n## Stages (by wall)\n\n");
    md.push_str("| Stage | Calls | Wall s | CPU s | Avg wall s | Alloc bytes |\n");
    md.push_str("|-------|------:|-------:|------:|-----------:|------------:|\n");
    let mut rows: Vec<_> = Stage::all()
        .iter()
        .map(|s| (*s, stages.get(*s)))
        .filter(|(_, st)| st.calls > 0 || st.alloc_bytes > 0)
        .collect();
    rows.sort_by(|a, b| b.1.wall_ns.cmp(&a.1.wall_ns));
    for (stage, st) in rows {
        md.push_str(&format!(
            "| `{}` | {} | {:.3} | {:.3} | {:.6} | {} |\n",
            stage.as_str(),
            st.calls,
            ns_to_s(st.wall_ns),
            ns_to_s(st.cpu_ns),
            st.avg_wall_ns() / 1e9,
            st.alloc_bytes
        ));
    }
    let occ = pairhmm.simd_occupancy();
    md.push_str("\n## PairHMM\n\n");
    md.push_str(&format!("- Reads scored: **{}**\n", pairhmm.reads_scored));
    md.push_str(&format!(
        "- Haplotypes scored: **{}**\n",
        pairhmm.haplotypes_scored
    ));
    md.push_str(&format!(
        "- Read×hap pairs: **{}**\n",
        pairhmm.read_hap_pairs
    ));
    md.push_str(&format!(
        "- Mean haps/read: **{:.2}**\n",
        pairhmm.mean_haps_per_read()
    ));
    md.push_str(&format!(
        "- SIMD pack units: **{}** ({}-wide ≈{} haps, occupancy **{:.1}%**)\n",
        occ.pack_units, occ.pack_lanes, occ.pack_hap_est, occ.pack_occupancy_pct
    ));
    md.push_str(&format!(
        "- Prefix reuse: **{:.1}%** ({} haps)\n",
        occ.prefix_reuse_pct, pairhmm.prefix_reuse_haps
    ));
    md.push_str(&format!(
        "- Leftover: **{:.1}%** ({} haps)\n",
        occ.leftover_pct, pairhmm.leftover_haps
    ));
    md.push_str(&format!(
        "- DP cells evaluated: **{}**; avoided via prefix: **{}**\n",
        pairhmm.dp_cells_evaluated, pairhmm.dp_cells_avoided_prefix
    ));
    md.push_str(&format!(
        "- PairHMM wall: **{:.3}s**\n",
        ns_to_s(pairhmm.wall_ns)
    ));

    md.push_str("\n## Genotyping\n\n");
    md.push_str(&format!("- Sites: **{}**\n", genotype.sites));
    md.push_str(&format!(
        "- Candidate alleles (sum): **{}**\n",
        genotype.candidate_alleles_sum
    ));
    md.push_str(&format!(
        "- Genotype states (sum): **{}** (diploid PL entries; typically 3/site)\n",
        genotype.genotype_states_sum
    ));
    md.push_str(&format!(
        "- PL vector lengths (sum): **{}**; samples (sum): **{}**\n",
        genotype.pl_vector_len_sum, genotype.samples_sum
    ));
    md.push_str(&format!(
        "- Time/site: **{:.6}s**; time/state: **{:.9}s**\n",
        genotype.time_per_site_ns() / 1e9,
        genotype.time_per_state_ns() / 1e9
    ));
    md.push_str(&format!(
        "- Nested — AD: **{:.3}s**; allele-map: **{:.3}s**; marginalize: **{:.3}s**; PL-enum: **{:.3}s**; event-rebuild: **{:.3}s**; assign Σ: **{:.3}s**\n",
        ns_to_s(genotype.ad_wall_ns),
        ns_to_s(genotype.allele_map_wall_ns),
        ns_to_s(genotype.marginalize_wall_ns),
        ns_to_s(genotype.genotype_enum_wall_ns),
        ns_to_s(genotype.event_rebuild_wall_ns),
        ns_to_s(genotype.wall_ns)
    ));
    md.push('\n');
    md
}

#[cfg(test)]
pub(super) fn build_json_for_test(
    stages: &StageAgg,
    pairhmm: &PairHmmAgg,
    genotype: &GenotypeAgg,
    run_wall: Duration,
    run_cpu: Option<Duration>,
    rayon_threads: usize,
    regions: u64,
) -> String {
    build_json(
        stages,
        pairhmm,
        genotype,
        run_wall,
        run_cpu,
        rayon_threads,
        regions,
    )
}
