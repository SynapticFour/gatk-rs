//! VCF record parity dumps (`j2-vcf`, `j2-format`).

use crate::engine::{CallRegionArgs, HaplotypeCallerEngine};
use crate::genotyping::emit_genotype_format_fields;
use crate::hc_genotyping_engine::{HcGenotypingConfig, RegionGenotypeResult};
use crate::read_model::ReadFilterParams;
use crate::region_vcf_emit::{
    build_biallelic_variant_record, try_emit_call_region_variant, try_emit_call_region_variants,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use crate::walker_apply::{call_disposition, AssemblyRegionCallDisposition};
use crate::walker_traversal::{
    flatten_assembly_regions, traverse_assembly_region_walker, WalkerTraversalConfig,
};
use gatk_common::{GatkError, GatkResult};
use gatk_core::io::vcf::{Genotype, VcfRecord};
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use std::io::Write;
use std::path::Path;

fn format_gt(gt: &Genotype) -> String {
    let a = gt.alleles.first().copied().unwrap_or(-1);
    let b = gt.alleles.get(1).copied().unwrap_or(-1);
    if a < 0 || b < 0 {
        "./.".to_string()
    } else {
        format!("{a}/{b}")
    }
}

pub fn write_variant_vcf_tsv(out: &mut impl Write, rec: Option<&VcfRecord>) -> GatkResult<()> {
    match rec {
        None => {
            writeln!(out, "variant_emitted\tfalse")?;
        }
        Some(r) => {
            writeln!(out, "variant_emitted\ttrue")?;
            writeln!(out, "chrom\t{}", r.chromosome)?;
            writeln!(out, "pos\t{}", r.position)?;
            writeln!(out, "id\t{}", r.id)?;
            writeln!(out, "ref\t{}", r.reference)?;
            writeln!(
                out,
                "alt\t{}",
                r.alternate
                    .iter()
                    .map(|a| a.as_str())
                    .collect::<Vec<_>>()
                    .join(",")
            )?;
            writeln!(
                out,
                "qual\t{}",
                r.quality
                    .map(|q| format!("{q:.6}"))
                    .unwrap_or_else(|| ".".to_string())
            )?;
            writeln!(
                out,
                "filter\t{}",
                r.filter
                    .iter()
                    .map(|s| s.as_str())
                    .collect::<Vec<_>>()
                    .join(";")
            )?;
        }
    }
    Ok(())
}

pub fn write_variant_format_tsv(out: &mut impl Write, rec: Option<&VcfRecord>) -> GatkResult<()> {
    let Some(r) = rec else {
        writeln!(out, "format_emitted\tfalse")?;
        return Ok(());
    };
    let Some(sample) = r.samples.first() else {
        writeln!(out, "format_emitted\tfalse")?;
        return Ok(());
    };
    writeln!(out, "format_emitted\ttrue")?;
    if let Some(gt) = &sample.gt {
        writeln!(out, "gt\t{}", format_gt(gt))?;
    }
    if let Some(gq) = sample.gq {
        writeln!(out, "gq\t{}", gq.round() as i32)?;
    }
    if let Some(dp) = sample.dp {
        writeln!(out, "dp\t{dp}")?;
    }
    if let Some(ad) = &sample.ad {
        let s = ad
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(",");
        writeln!(out, "ad\t{s}")?;
    }
    if let Some(pl) = &sample.pl {
        let s = pl
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(",");
        writeln!(out, "pl\t{s}")?;
    }
    Ok(())
}

/// First active region `call_region` → VCF identity fields (`j2-vcf` / `call-region-vcf`).
pub fn dump_call_region_vcf_tsv(
    ref_fasta: &Path,
    bam_path: &Path,
    interval_cli: &str,
    padding: u64,
    out: &mut impl Write,
) -> GatkResult<()> {
    let dict = SequenceDictionary::from_fasta_path(ref_fasta)?;
    let specs = parse_intervals_cli_string(&dict, interval_cli)?;
    let filters = ReadFilterParams::gatk_standard_hc();
    let cfg = WalkerTraversalConfig::gatk_haplotype_caller_production(padding);
    let walk = traverse_assembly_region_walker(&dict, &specs, ref_fasta, bam_path, &filters, &cfg)?;
    let regions = flatten_assembly_regions(&walk);
    let args = CallRegionArgs::default();
    let mut emitted: Option<VcfRecord> = None;
    for region in regions.iter().filter(|r| {
        matches!(
            call_disposition(r),
            AssemblyRegionCallDisposition::ActiveFull
        )
    }) {
        let Some(outcome) = HaplotypeCallerEngine::call_region(region, &dict, ref_fasta, &args)?
        else {
            continue;
        };
        if let Some(rec) = try_emit_call_region_variants(
            region,
            &outcome,
            "SAMPLE",
            DEFAULT_STAND_EMIT_CONFIDENCE,
        )?
        .into_iter()
        .next()
        {
            emitted = Some(rec);
            break;
        }
    }
    write_variant_vcf_tsv(out, emitted.as_ref())
}

/// Deterministic SNP from p7-style GL/AD fixture (`j2-vcf` synthetic row).
pub fn dump_variant_vcf_from_gl_ad_tsv(
    contig: &str,
    pos_1based: u64,
    ref_allele: &str,
    alt_allele: &str,
    gl_csv: &str,
    ad_csv: &str,
    out: &mut impl Write,
) -> GatkResult<()> {
    let gls: Vec<f64> = gl_csv
        .split(',')
        .map(|s| s.parse())
        .collect::<Result<_, _>>()
        .map_err(|e| GatkError::argument(format!("parse gl: {e}")))?;
    let ads: Vec<i32> = ad_csv
        .split(',')
        .map(|s| s.parse())
        .collect::<Result<_, _>>()
        .map_err(|e| GatkError::argument(format!("parse ad: {e}")))?;
    let format = emit_genotype_format_fields(&gls, &ads)?;
    let genotype = RegionGenotypeResult {
        aggregation: crate::genotyping::HaplotypeLikelihoodAggregation {
            haplotype_log10_sums: vec![0.0, -1.0],
            read_count: 1,
        },
        best_haplotype_index: 0,
        ref_haplotype_index: 0,
        alt_haplotype_index: 1,
        genotype_log10_likelihoods: gls,
        format,
    };
    let rec = build_biallelic_variant_record(
        contig,
        None,
        pos_1based,
        ref_allele,
        alt_allele,
        &genotype,
        "SAMPLE",
        &HcGenotypingConfig::default(),
    )?;
    let gt = rec
        .samples
        .first()
        .and_then(|s| s.gt.as_ref())
        .map(|g| g.alleles.as_slice())
        .unwrap_or(&[]);
    let emit = !matches!(gt, [0, 0])
        && genotype.format.gq.as_i32() as f64 >= DEFAULT_STAND_EMIT_CONFIDENCE;
    write_variant_vcf_tsv(out, if emit { Some(&rec) } else { None })
}

pub fn dump_call_region_format_tsv(
    ref_fasta: &Path,
    bam_path: &Path,
    interval_cli: &str,
    padding: u64,
    out: &mut impl Write,
) -> GatkResult<()> {
    let dict = SequenceDictionary::from_fasta_path(ref_fasta)?;
    let specs = parse_intervals_cli_string(&dict, interval_cli)?;
    let filters = ReadFilterParams::gatk_standard_hc();
    let cfg = WalkerTraversalConfig::gatk_haplotype_caller_production(padding);
    let walk = traverse_assembly_region_walker(&dict, &specs, ref_fasta, bam_path, &filters, &cfg)?;
    let regions = flatten_assembly_regions(&walk);
    let args = CallRegionArgs::default();
    let mut emitted: Option<VcfRecord> = None;
    for region in regions.iter().filter(|r| {
        matches!(
            call_disposition(r),
            AssemblyRegionCallDisposition::ActiveFull
        )
    }) {
        let Some(outcome) = HaplotypeCallerEngine::call_region(region, &dict, ref_fasta, &args)?
        else {
            continue;
        };
        if let Some(rec) =
            try_emit_call_region_variant(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)?
        {
            emitted = Some(rec);
            break;
        }
    }
    write_variant_format_tsv(out, emitted.as_ref())
}

pub fn dump_variant_format_from_gl_ad_tsv(
    contig: &str,
    pos_1based: u64,
    ref_allele: &str,
    alt_allele: &str,
    gl_csv: &str,
    ad_csv: &str,
    out: &mut impl Write,
) -> GatkResult<()> {
    let mut buf = Vec::new();
    dump_variant_vcf_from_gl_ad_tsv(
        contig, pos_1based, ref_allele, alt_allele, gl_csv, ad_csv, &mut buf,
    )?;
    let gls: Vec<f64> = gl_csv
        .split(',')
        .map(|s| s.parse())
        .collect::<Result<_, _>>()
        .map_err(|e| GatkError::argument(format!("parse gl: {e}")))?;
    let ads: Vec<i32> = ad_csv
        .split(',')
        .map(|s| s.parse())
        .collect::<Result<_, _>>()
        .map_err(|e| GatkError::argument(format!("parse ad: {e}")))?;
    let format = emit_genotype_format_fields(&gls, &ads)?;
    let genotype = RegionGenotypeResult {
        aggregation: crate::genotyping::HaplotypeLikelihoodAggregation {
            haplotype_log10_sums: vec![0.0, -1.0],
            read_count: 1,
        },
        best_haplotype_index: 0,
        ref_haplotype_index: 0,
        alt_haplotype_index: 1,
        genotype_log10_likelihoods: gls,
        format,
    };
    let rec = build_biallelic_variant_record(
        contig,
        None,
        pos_1based,
        ref_allele,
        alt_allele,
        &genotype,
        "SAMPLE",
        &HcGenotypingConfig::default(),
    )?;
    let gt = rec
        .samples
        .first()
        .and_then(|s| s.gt.as_ref())
        .map(|g| g.alleles.as_slice())
        .unwrap_or(&[]);
    let emit = !matches!(gt, [0, 0])
        && genotype.format.gq.as_i32() as f64 >= DEFAULT_STAND_EMIT_CONFIDENCE;
    write_variant_format_tsv(out, if emit { Some(&rec) } else { None })
}
