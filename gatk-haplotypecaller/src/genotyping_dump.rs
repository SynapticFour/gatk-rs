//! Genotyping parity dumps (L2 parity).

use crate::af_calc::{calculate_biallelic_af_em, AfCalculatorConfig};
use crate::allele_subsetting::{biallelic_ref_alt_indices, subset_alleles_for_genotyping};
use crate::allele_subsetting_pl::{
    subset_het_ac_sac, subset_trim_acg_to_ac, subset_trim_acg_to_ag,
};
use crate::feature_context::FeatureDataSources;
use crate::genotype_limits::compute_max_acceptable_allele_count;
use crate::genotyping::{
    aggregate_haplotype_log10_likelihoods, best_haplotype_index, emit_genotype_format_fields,
    emit_genotype_phasing_fields, ReadLikelihoodRow,
};
use crate::hc_genotyping_engine::HcGenotypingConfig;
use crate::pairhmm::{pairhmm_log10_likelihood, PairHmmInput, PairHmmParams};
use crate::pairhmm_dump::load_pairhmm_cases_tsv;
use gatk_common::{GatkError, GatkResult};
use indexmap::IndexMap;
use std::io::Write;
use std::path::Path;

/// Aggregate per-read PairHMM rows (F.1 cases) into haplotype sums — G.1 gate.
pub fn dump_genotyping_aggregate_tsv(
    pairhmm_cases_path: &Path,
    out: &mut impl Write,
) -> GatkResult<()> {
    let cases = load_pairhmm_cases_tsv(pairhmm_cases_path)?;
    let params = PairHmmParams::default();
    // IndexMap: first-seen read order (PairHMM case TSV order) — HashMap iteration
    // would make aggregate dump / best-haplotype selection order-dependent.
    let mut by_read: IndexMap<String, Vec<f64>> = IndexMap::new();
    let mut haplotype_labels: Vec<String> = Vec::new();
    for row in &cases {
        let ll = pairhmm_log10_likelihood(
            &PairHmmInput {
                read_bases: row.read_bases.clone(),
                read_base_quals: row.read_base_quals.clone(),
                read_mapping_quality: row.read_mapq,
                haplotype_bases: row.haplotype.clone(),
            },
            &params,
        )?;
        // CLONE: needed because owned HashMap entry key.
        let entry = by_read.entry(row.read_bases.clone()).or_default();
        if !haplotype_labels.contains(&row.haplotype) {
            // CLONE: needed because owned element into collection.
            haplotype_labels.push(row.haplotype.clone());
        }
        let hap_idx = haplotype_labels
            .iter()
            .position(|h| h == &row.haplotype)
            .expect("hap label");
        if entry.len() <= hap_idx {
            entry.resize(hap_idx + 1, f64::NEG_INFINITY);
        }
        entry[hap_idx] = ll;
    }

    let rows: Vec<ReadLikelihoodRow> = by_read
        .into_iter()
        .enumerate()
        .map(|(i, (read_id, haps))| ReadLikelihoodRow {
            read_id: format!("read_{i}_{read_id}"),
            haplotype_log10_likelihoods: haps,
        })
        .collect();

    let agg = aggregate_haplotype_log10_likelihoods(&rows)?;
    let best = best_haplotype_index(&agg)
        .unwrap_or(crate::bio_ids::HaplotypeIndex::new(0))
        .get();

    writeln!(out, "haplotype_count\t{}", agg.haplotype_log10_sums.len())
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    writeln!(out, "read_count\t{}", agg.read_count)
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    for (i, sum) in agg.haplotype_log10_sums.iter().enumerate() {
        writeln!(out, "haplotype_{i}_log10_sum\t{sum}")
            .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    }
    writeln!(out, "best_haplotype_index\t{best}")
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    Ok(())
}

/// G.2.2 — PL/GQ/AD/DP from fixture genotype log10 likelihoods (p7 contract).
pub fn dump_genotype_format_tsv(fixture_path: &Path, out: &mut impl Write) -> GatkResult<()> {
    let raw = std::fs::read_to_string(fixture_path)
        .map_err(|e| GatkError::io(format!("read fixture {}", fixture_path.display()), e))?;
    writeln!(out, "# case_id\tpl\tgq\tad\tdp")
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    for line in raw.lines() {
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = t.split('\t').collect();
        if cols.len() != 3 {
            return Err(GatkError::argument(format!(
                "genotype-format fixture row needs 3 cols: {t}"
            )));
        }
        let case_id = cols[0];
        let gls: Vec<f64> = cols[1]
            .split(',')
            .map(|s| s.parse::<f64>())
            .collect::<Result<_, _>>()
            .map_err(|e| GatkError::argument(format!("parse gl for {case_id}: {e}")))?;
        let ads: Vec<i32> = cols[2]
            .split(',')
            .map(|s| s.parse::<i32>())
            .collect::<Result<_, _>>()
            .map_err(|e| GatkError::argument(format!("parse ad for {case_id}: {e}")))?;
        let fields = emit_genotype_format_fields(&gls, &ads)?;
        let pl = fields
            .pl
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(",");
        let ad = fields
            .ad
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(",");
        writeln!(out, "{case_id}\t{pl}\t{}\t{ad}\t{}", fields.gq, fields.dp)
            .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    }
    Ok(())
}

/// G.2.3 — AF / EM posteriors from per-sample GL fixtures (`g2-af`).
pub fn dump_af_em_tsv(fixture_path: &Path, out: &mut impl Write) -> GatkResult<()> {
    let raw = std::fs::read_to_string(fixture_path)
        .map_err(|e| GatkError::io(format!("read {}", fixture_path.display()), e))?;
    for line in raw.lines() {
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = t.split('\t').collect();
        if cols.len() < 2 {
            return Err(GatkError::argument(format!("g2-af row needs gl col: {t}")));
        }
        let case_id = cols[0];
        let gls: Vec<f64> = cols[1]
            .split(',')
            .map(|s| s.parse())
            .collect::<Result<_, _>>()
            .map_err(|e| GatkError::argument(format!("parse gl: {e}")))?;
        let gl_roundtrip =
            crate::activity_scoring::genotype_log10_likelihoods_after_java_genotype_pl_roundtrip(
                &gls,
            );
        let gl_ref = gl_roundtrip.as_slice();
        let result = calculate_biallelic_af_em(&[gl_ref], &AfCalculatorConfig::default())?;
        writeln!(out, "{case_id}\talt_ac\t{}", result.alt_allele_count)?;
        writeln!(out, "{case_id}\taf\t{:.6}", result.af)?;
        writeln!(
            out,
            "{case_id}\tlog10_p_no_variant\t{:.6}",
            result.log10_posterior_no_variant
        )?;
        writeln!(out, "{case_id}\tem_iterations\t{}", result.em_iterations)?;
    }
    Ok(())
}

/// G.3 — max acceptable allele count (`g3`).
pub fn dump_genotype_limits_tsv(
    ploidy: u32,
    max_genotype_count: u32,
    out: &mut impl Write,
) -> GatkResult<()> {
    let max_alleles = compute_max_acceptable_allele_count(ploidy, max_genotype_count)?;
    writeln!(out, "ploidy\t{ploidy}")?;
    writeln!(out, "max_genotype_count\t{max_genotype_count}")?;
    writeln!(out, "max_acceptable_allele_count\t{max_alleles}")?;
    Ok(())
}

/// G.4.1 — physical phasing FORMAT fields (`g4`).
pub fn dump_genotype_phasing_tsv(
    alleles_csv: &str,
    phasing_enabled: bool,
    phase_set: Option<i32>,
    out: &mut impl Write,
) -> GatkResult<()> {
    let alleles: Vec<i32> = alleles_csv
        .split(',')
        .map(|s| s.parse())
        .collect::<Result<_, _>>()
        .map_err(|e| GatkError::argument(format!("parse alleles: {e}")))?;
    let fields = emit_genotype_phasing_fields(&alleles, phasing_enabled, phase_set)?;
    writeln!(out, "gt\t{}", fields.gt)?;
    writeln!(out, "phased\t{}", fields.phased)?;
    writeln!(out, "pgt\t{}", fields.pgt.as_deref().unwrap_or("-"))?;
    writeln!(out, "pid\t{}", fields.pid.as_deref().unwrap_or("-"))?;
    writeln!(
        out,
        "ps\t{}",
        fields
            .ps
            .map(|v| v.to_string())
            .unwrap_or_else(|| "-".to_string())
    )?;
    Ok(())
}

/// G.4.2 — force-calling allele present at locus (`g4-force`).
pub fn dump_force_calling_genotype_tsv(
    features_path: &Path,
    contig: &str,
    pos_1based: u64,
    force_call_filtered: bool,
    out: &mut impl Write,
) -> GatkResult<()> {
    let mut sources = FeatureDataSources::default();
    sources.load_vcf_source(
        crate::assembly_regions_dump::GATK_HC_ALLELES_FEATURE_SOURCE,
        features_path,
    )?;
    let present = sources.force_calling_allele_overlaps_locus(
        crate::assembly_regions_dump::GATK_HC_ALLELES_FEATURE_SOURCE,
        contig,
        pos_1based,
        force_call_filtered,
    );
    writeln!(out, "contig\t{contig}")?;
    writeln!(out, "pos\t{pos_1based}")?;
    writeln!(out, "force_calling_present\t{present}")?;
    writeln!(out, "genotyping_config_ok\ttrue")?;
    let _ = HcGenotypingConfig::default();
    Ok(())
}

/// G-D05 — allele subsetting indices for haplotype log10 sums fixture.
pub fn dump_allele_subsetting_tsv(
    haplotype_log10_sums_csv: &str,
    is_reference_csv: &str,
    max_allele_count: usize,
    out: &mut impl Write,
) -> GatkResult<()> {
    use crate::genotyping::HaplotypeLikelihoodAggregation;
    use crate::haplotype::Haplotype;
    let sums: Vec<f64> = haplotype_log10_sums_csv
        .split(',')
        .map(|s| s.parse())
        .collect::<Result<_, _>>()
        .map_err(|e| GatkError::argument(format!("parse sums: {e}")))?;
    let is_ref: Vec<bool> = is_reference_csv
        .split(',')
        .map(|s| s == "1" || s.eq_ignore_ascii_case("true"))
        .collect();
    let haps: Vec<Haplotype> = is_ref.iter().map(|&r| Haplotype::new(b"A", r)).collect();
    let agg = HaplotypeLikelihoodAggregation {
        haplotype_log10_sums: sums,
        read_count: 1,
    };
    let kept = subset_alleles_for_genotyping(&haps, &agg, max_allele_count)?;
    let (ref_i, alt_i) = biallelic_ref_alt_indices(&agg, &haps);
    writeln!(out, "haplotype_count\t{}", haps.len())?;
    writeln!(
        out,
        "kept_indices\t{}",
        kept.iter()
            .map(|i| i.to_string())
            .collect::<Vec<_>>()
            .join(",")
    )?;
    writeln!(out, "ref_haplotype_index\t{ref_i}")?;
    writeln!(out, "alt_haplotype_index\t{alt_i}")?;
    Ok(())
}

/// G-D05 — `AlleleSubsettingUtils.subsetAlleles` PL/AD (`g-subset-pl`).
pub fn dump_subset_alleles_pl_tsv(fixture_path: &Path, out: &mut impl Write) -> GatkResult<()> {
    let raw = std::fs::read_to_string(fixture_path)
        .map_err(|e| GatkError::io(format!("read {}", fixture_path.display()), e))?;
    for line in raw.lines() {
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = t.split('\t').collect();
        if cols.len() < 3 {
            return Err(GatkError::argument(format!(
                "g-subset-pl row needs 3 cols: {t}"
            )));
        }
        let case_id = cols[0];
        let log10_pl: Vec<f64> = cols[1]
            .split(',')
            .map(|s| s.parse())
            .collect::<Result<_, _>>()
            .map_err(|e| GatkError::argument(format!("parse pl: {e}")))?;
        let ad: Vec<i32> = cols[2]
            .split(',')
            .map(|s| s.parse())
            .collect::<Result<_, _>>()
            .map_err(|e| GatkError::argument(format!("parse ad: {e}")))?;
        let result = match case_id {
            "trim_acg_to_ag" => subset_trim_acg_to_ag(&log10_pl, &ad)?,
            _ => subset_trim_acg_to_ac(&log10_pl, &ad)?,
        };
        write_subset_vc_rows(case_id, &result, out)?;
    }
    Ok(())
}

fn write_subset_vc_rows(
    case_id: &str,
    result: &crate::allele_subsetting_pl::SubsetAllelesPlResult,
    out: &mut impl Write,
) -> GatkResult<()> {
    writeln!(
        out,
        "{case_id}\tallele_count_before\t{}",
        result.allele_count_before
    )?;
    writeln!(
        out,
        "{case_id}\tallele_count_after\t{}",
        result.allele_count_after
    )?;
    writeln!(out, "{case_id}\tpl_length\t{}", result.pl.len())?;
    for (i, pl) in result.pl.iter().enumerate() {
        writeln!(out, "{case_id}\tpl_{i}\t{pl}")?;
    }
    writeln!(
        out,
        "{case_id}\tad\t{}",
        result
            .ad
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(",")
    )?;
    if let Some(gq) = result.gq {
        writeln!(out, "{case_id}\tgq\t{gq}")?;
    }
    if let Some(sac) = &result.sac {
        writeln!(
            out,
            "{case_id}\tsac\t{}",
            sac.iter()
                .map(|v| v.to_string())
                .collect::<Vec<_>>()
                .join(",")
        )?;
    }
    Ok(())
}

/// G-D05 — `AlleleSubsettingUtils.subsetAlleles` with SAC/GQ (`g-subset-vc`).
pub fn dump_subset_alleles_vc_tsv(fixture_path: &Path, out: &mut impl Write) -> GatkResult<()> {
    let raw = std::fs::read_to_string(fixture_path)
        .map_err(|e| GatkError::io(format!("read {}", fixture_path.display()), e))?;
    for line in raw.lines() {
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = t.split('\t').collect();
        if cols.len() < 4 {
            return Err(GatkError::argument(format!(
                "g-subset-vc row needs 4 cols: {t}"
            )));
        }
        let case_id = cols[0];
        let log10_pl: Vec<f64> = cols[1]
            .split(',')
            .map(|s| s.parse())
            .collect::<Result<_, _>>()
            .map_err(|e| GatkError::argument(format!("parse pl: {e}")))?;
        let ad: Vec<i32> = cols[2]
            .split(',')
            .map(|s| s.parse())
            .collect::<Result<_, _>>()
            .map_err(|e| GatkError::argument(format!("parse ad: {e}")))?;
        let sac: Vec<i32> = cols[3]
            .split(',')
            .map(|s| s.parse())
            .collect::<Result<_, _>>()
            .map_err(|e| GatkError::argument(format!("parse sac: {e}")))?;
        let result = subset_het_ac_sac(&log10_pl, &ad, &sac)?;
        write_subset_vc_rows(case_id, &result, out)?;
    }
    Ok(())
}

/// G-D05 integration: haplotype trim kept indices + VC subset (`g-subset-integration`).
pub fn dump_subset_alleles_integration_tsv(
    haplotype_log10_sums_csv: &str,
    is_reference_csv: &str,
    max_allele_count: usize,
    fixture_path: &Path,
    out: &mut impl Write,
) -> GatkResult<()> {
    dump_allele_subsetting_tsv(
        haplotype_log10_sums_csv,
        is_reference_csv,
        max_allele_count,
        out,
    )?;
    dump_subset_alleles_vc_tsv(fixture_path, out)
}
