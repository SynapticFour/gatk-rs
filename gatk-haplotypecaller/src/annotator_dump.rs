//! Annotated variant site parity dumps.

use crate::annotator::{annotate_parity_v1_site, AnnotatedSite, VariantAnnotationContext};
use crate::genotyping::SampleAnnotationInput;
use gatk_common::{GatkError, GatkResult};
use std::io::Write;
use std::path::Path;

fn parse_samples_tsv(path: &Path) -> GatkResult<Vec<SampleAnnotationInput>> {
    let raw = std::fs::read_to_string(path)
        .map_err(|e| GatkError::io(format!("read {}", path.display()), e))?;
    let mut out = Vec::new();
    for line in raw.lines() {
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = t.split('\t').collect();
        if cols.len() < 2 {
            return Err(GatkError::argument(format!(
                "samples row needs >=2 cols: {t}"
            )));
        }
        let gts: Vec<i32> = cols[1]
            .split(',')
            .map(|s| s.parse::<i32>())
            .collect::<Result<_, _>>()
            .map_err(|e| GatkError::argument(format!("parse genotype alleles: {e}")))?;
        let dp = if cols.len() >= 3 && !cols[2].is_empty() && cols[2] != "-" {
            Some(
                cols[2]
                    .parse::<i32>()
                    .map_err(|e| GatkError::argument(format!("parse dp: {e}")))?,
            )
        } else {
            None
        };
        out.push(SampleAnnotationInput {
            genotype_alleles: gts,
            dp,
        });
    }
    Ok(out)
}

fn write_annotated_site(out: &mut impl Write, site: &AnnotatedSite) -> GatkResult<()> {
    writeln!(out, "info_key_count\t{}", site.info_keys.len())?;
    for (i, k) in site.info_keys.iter().enumerate() {
        writeln!(out, "info_key_{i}\t{k}")?;
    }
    writeln!(out, "format_key_count\t{}", site.format_keys.len())?;
    for (i, k) in site.format_keys.iter().enumerate() {
        writeln!(out, "format_key_{i}\t{k}")?;
    }
    writeln!(out, "ac\t{}", format_int_list(&site.ac))?;
    writeln!(out, "an\t{}", site.an)?;
    writeln!(out, "af\t{}", format_float_list(&site.af))?;
    writeln!(out, "ns\t{}", site.ns)?;
    writeln!(out, "dp\t{}", site.dp)?;
    Ok(())
}

fn format_int_list(v: &[i32]) -> String {
    v.iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn format_float_list(v: &[f64]) -> String {
    v.iter()
        .map(|x| format!("{x:.6}"))
        .collect::<Vec<_>>()
        .join(",")
}

/// Run parity v1 annotator on a samples fixture and emit stable TSV (`i1-core` gate).
pub fn dump_annotate_core_tsv(
    alt_allele_count: usize,
    samples_path: &Path,
    out: &mut impl Write,
) -> GatkResult<()> {
    let samples = parse_samples_tsv(samples_path)?;
    let ctx = VariantAnnotationContext {
        alt_allele_count,
        samples,
    };
    let site = annotate_parity_v1_site(&ctx)?;
    write_annotated_site(out, &site)
}

/// I-D01 — standard annotation plugin slice dump.
pub fn dump_standard_annotations_tsv(
    ref_fw: u32,
    ref_rv: u32,
    alt_fw: u32,
    alt_rv: u32,
    qual: f64,
    dp: i32,
    ref_bqs_csv: &str,
    alt_bqs_csv: &str,
    ref_pos_csv: &str,
    alt_pos_csv: &str,
    ref_mq_csv: &str,
    alt_mq_csv: &str,
    out: &mut impl Write,
) -> GatkResult<()> {
    use crate::annotator::plugins::{
        fisher_strand, mapping_quality_rank_sum, qual_by_depth, rank_sum_baseq, read_pos_rank_sum,
        strand_odds_ratio,
    };
    let parse_u8 = |s: &str| -> Vec<u8> {
        if s.is_empty() || s == "-" {
            return Vec::new();
        }
        s.split(',').filter_map(|x| x.trim().parse().ok()).collect()
    };
    let parse_f64 = |s: &str| -> Vec<f64> {
        if s.is_empty() || s == "-" {
            return Vec::new();
        }
        s.split(',').filter_map(|x| x.trim().parse().ok()).collect()
    };
    let fs = fisher_strand::fisher_strand_statistic(ref_fw, ref_rv, alt_fw, alt_rv);
    let sor = strand_odds_ratio::strand_odds_ratio(ref_fw, ref_rv, alt_fw, alt_rv);
    let qd = qual_by_depth::qual_by_depth(qual, dp);
    let bq = rank_sum_baseq::base_quality_rank_sum(&parse_u8(ref_bqs_csv), &parse_u8(alt_bqs_csv));
    let rp = read_pos_rank_sum::read_pos_rank_sum(&parse_f64(ref_pos_csv), &parse_f64(alt_pos_csv));
    let mq = mapping_quality_rank_sum::mapping_quality_rank_sum(
        &parse_u8(ref_mq_csv),
        &parse_u8(alt_mq_csv),
    );
    writeln!(out, "FS\t{fs:.6}")?;
    writeln!(out, "SOR\t{sor:.6}")?;
    writeln!(out, "QD\t{qd:.6}")?;
    writeln!(out, "BaseQRankSum\t{bq:.6}")?;
    writeln!(out, "ReadPosRankSum\t{rp:.6}")?;
    writeln!(out, "MQRankSum\t{mq:.6}")?;
    Ok(())
}

/// I-D02 — allele-specific INFO scaffold.
pub fn dump_as_annotations_tsv(
    site_af: f64,
    site_qual: f64,
    out: &mut impl Write,
) -> GatkResult<()> {
    use crate::annotator::plugins::as_standard;
    writeln!(out, "AS_AF_0\t{:.6}", as_standard::as_af(site_af, 0))?;
    writeln!(out, "AS_AF_1\t{:.6}", as_standard::as_af(site_af, 1))?;
    writeln!(
        out,
        "AS_QUAL_0\t{:.6}",
        as_standard::as_qual(site_qual, site_af, 0)
    )?;
    writeln!(
        out,
        "AS_QUAL_1\t{:.6}",
        as_standard::as_qual(site_qual, site_af, 1)
    )?;
    Ok(())
}

/// I-D03 — ExcessHet (GATK `ExcessHet.calculateEH` on genotype counts).
pub fn dump_excess_het_tsv(
    ref_count: u32,
    het_count: u32,
    hom_count: u32,
    out: &mut impl Write,
) -> GatkResult<()> {
    use crate::annotator::plugins::excess_het;
    let eh = excess_het::excess_heterozygosity_phred(ref_count, het_count, hom_count);
    writeln!(out, "ExcessHet\t{eh:.6}")?;
    Ok(())
}

/// I-D04 — DepthPerSampleHC vs FORMAT DP.
pub fn dump_depth_per_sample_hc_tsv(ad_csv: &str, out: &mut impl Write) -> GatkResult<()> {
    use crate::annotator::plugins::depth_per_sample_hc;
    let ad: Vec<i32> = if ad_csv.is_empty() || ad_csv == "-" {
        Vec::new()
    } else {
        ad_csv
            .split(',')
            .map(|s| s.trim().parse())
            .collect::<Result<_, _>>()
            .map_err(|e| GatkError::argument(format!("parse ad: {e}")))?
    };
    let dps = depth_per_sample_hc::depth_per_sample_hc(&ad);
    let dp = depth_per_sample_hc::format_dp_from_ad(&ad);
    writeln!(out, "DepthPerSampleHC\t{dps}")?;
    writeln!(out, "FORMAT_DP\t{dp}")?;
    writeln!(out, "reconciled\t{}", dps == dp)?;
    Ok(())
}

/// I-D05 — per-plugin gate rows (FS / QD / BaseQRankSum).
pub fn dump_annotation_plugin_tsv(
    plugin: &str,
    ref_fw: u32,
    ref_rv: u32,
    alt_fw: u32,
    alt_rv: u32,
    qual: f64,
    dp: i32,
    ref_bqs_csv: &str,
    alt_bqs_csv: &str,
    out: &mut impl Write,
) -> GatkResult<()> {
    use crate::annotator::plugins::{fisher_strand, qual_by_depth, rank_sum_baseq};
    let parse_bqs = |s: &str| -> Vec<u8> {
        if s.is_empty() || s == "-" {
            return Vec::new();
        }
        s.split(',').filter_map(|x| x.trim().parse().ok()).collect()
    };
    match plugin {
        "FS" | "fisher_strand" => {
            let v = fisher_strand::fisher_strand_statistic(ref_fw, ref_rv, alt_fw, alt_rv);
            writeln!(out, "plugin\tFS")?;
            writeln!(out, "value\t{v:.6}")?;
        }
        "QD" | "qual_by_depth" => {
            let v = qual_by_depth::qual_by_depth(qual, dp);
            writeln!(out, "plugin\tQD")?;
            writeln!(out, "value\t{v:.6}")?;
        }
        "BaseQRankSum" | "rank_sum_baseq" => {
            let v = rank_sum_baseq::base_quality_rank_sum(
                &parse_bqs(ref_bqs_csv),
                &parse_bqs(alt_bqs_csv),
            );
            writeln!(out, "plugin\tBaseQRankSum")?;
            writeln!(out, "value\t{v:.6}")?;
        }
        other => {
            return Err(GatkError::argument(format!(
                "unknown annotation plugin: {other}"
            )));
        }
    }
    Ok(())
}

/// Dump manifest line counts for parity v1 (I.0.1 sanity).
pub fn dump_annotation_manifest_tsv(out: &mut impl Write) -> GatkResult<()> {
    use crate::annotator::{PARITY_V1_FORMAT_KEYS, PARITY_V1_INFO_KEYS};
    writeln!(out, "parity_v1_info_count\t{}", PARITY_V1_INFO_KEYS.len())?;
    for (i, k) in PARITY_V1_INFO_KEYS.iter().enumerate() {
        writeln!(out, "parity_v1_info_{i}\t{k}")?;
    }
    writeln!(
        out,
        "parity_v1_format_count\t{}",
        PARITY_V1_FORMAT_KEYS.len()
    )?;
    for (i, k) in PARITY_V1_FORMAT_KEYS.iter().enumerate() {
        writeln!(out, "parity_v1_format_{i}\t{k}")?;
    }
    Ok(())
}
