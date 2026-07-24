//! Java `HcParityRegionGenotype` + `HcParityNativePairHmm` path for L2 parity G.2 live gates.
//! Uses raw shard reads (sorted by qname, then start), parity PairHMM defaults (45/45/10),
//! no finalize/overlap correction and no MQ/BQ capping — matching `HcFullParityGateDump.assemblyRegionGenotypeSubset`.

use crate::genotyping::{
    aggregate_haplotype_log10_likelihoods, emit_genotype_format_fields, GenotypeFormatFields,
    ReadLikelihoodRow,
};
use crate::haplotype::Haplotype;
use crate::hc_genotyping_engine::{
    biallelic_genotype_log10_likelihoods_gatk, biallelic_genotype_log10_likelihoods_parity_legacy,
};
use crate::pairhmm_log10::log10_pairhmm_likelihood_parity_defaults;
use gatk_common::{GatkError, GatkResult};
use std::io::Write;

/// Build per-read haplotype log10 rows like Java `assemblyRegionGenotypeSubset`.
pub fn parity_java_aligned_read_rows(
    region_reads: &[rust_htslib::bam::Record],
    haplotypes: &[Haplotype],
) -> GatkResult<Vec<ReadLikelihoodRow>> {
    let mut records: Vec<&rust_htslib::bam::Record> = region_reads.iter().collect();
    records.sort_by(|a, b| a.qname().cmp(b.qname()).then_with(|| a.pos().cmp(&b.pos())));

    let n_haps = haplotypes.len();
    let mut rows = Vec::with_capacity(records.len());
    for rec in records {
        let bases = rec.seq().as_bytes();
        let quals = rec.qual();
        let mut ll = vec![f64::NEG_INFINITY; n_haps];
        for (hi, hap) in haplotypes.iter().enumerate() {
            ll[hi] = log10_pairhmm_likelihood_parity_defaults(&bases, quals, hap.bases.as_slice())?;
        }
        let read_id = String::from_utf8_lossy(rec.qname()).into_owned();
        rows.push(ReadLikelihoodRow {
            read_id,
            haplotype_log10_likelihoods: ll,
        });
    }
    Ok(rows)
}

/// Sum read×hap matrix columns (Java live-subset hap score input).
pub fn parity_java_aligned_hap_log10_sums(
    rows: &[ReadLikelihoodRow],
) -> GatkResult<crate::genotyping::HaplotypeLikelihoodAggregation> {
    aggregate_haplotype_log10_likelihoods(rows)
}

/// Java `HcParityRegionGenotype.GenotypeDump` from per-read likelihood rows.
/// # Invariants
/// When `genotyped` is false, `genotype_log10` is empty and `format` is `None`.
/// `ref_hap_idx` / `alt_hap_idx` index into the haplotype list used for biallelic GLs.
/// # Ownership
/// Owns genotype vectors and optional [`GenotypeFormatFields`]; rows are consumed upstream.
/// # Mutation
/// Immutable snapshot for parity dumps after genotype computation.
/// # Biological assumptions
/// Biallelic ref/alt genotype from summed read×haplogtype log10 likelihoods (parity gate path).
/// # Java equivalence
/// Java `HcParityRegionGenotype.GenotypeDump` / `assemblyRegionGenotypeSubset` live gates.
pub struct ParityRegionGenotypeDump {
    pub haplotype_count: usize,
    pub read_count: usize,
    pub genotyped: bool,
    pub ref_hap_idx: usize,
    pub alt_hap_idx: usize,
    pub best_hap_idx: usize,
    pub genotype_log10: Vec<f64>,
    pub format: Option<GenotypeFormatFields>,
}

pub fn parity_region_genotype_from_rows(
    rows: &[ReadLikelihoodRow],
    is_reference: &[bool],
) -> GatkResult<ParityRegionGenotypeDump> {
    parity_region_genotype_from_rows_with_gl_mode(rows, is_reference, false)
}

/// `legacy_gl`: frozen `g2-region` Java dumps use `2×lr` / `lr+la` / `2×la` (no ploidy denominator).
pub fn parity_region_genotype_from_rows_with_gl_mode(
    rows: &[ReadLikelihoodRow],
    is_reference: &[bool],
    legacy_gl: bool,
) -> GatkResult<ParityRegionGenotypeDump> {
    let hap_count = is_reference.len();
    let read_count = rows.len();
    if hap_count == 0 || read_count == 0 {
        return Ok(ParityRegionGenotypeDump {
            haplotype_count: hap_count,
            read_count,
            genotyped: false,
            ref_hap_idx: 0,
            alt_hap_idx: 0,
            best_hap_idx: 0,
            genotype_log10: Vec::new(),
            format: None,
        });
    }
    let mut sums = vec![0.0_f64; hap_count];
    for row in rows {
        for (i, s) in sums.iter_mut().enumerate() {
            *s += row.haplotype_log10_likelihoods[i];
        }
    }
    let mut best = 0usize;
    for i in 1..hap_count {
        if sums[i] > sums[best] {
            best = i;
        }
    }
    let ref_idx = is_reference.iter().position(|&r| r).unwrap_or(0);
    let mut alt_idx = ref_idx;
    let mut alt_sum = f64::NEG_INFINITY;
    for (i, &is_ref) in is_reference.iter().enumerate() {
        if i != ref_idx && !is_ref && sums[i] > alt_sum {
            alt_sum = sums[i];
            alt_idx = i;
        }
    }
    let mut ref_ad = 0i32;
    let mut alt_ad = 0i32;
    for row in rows {
        let lr = row.haplotype_log10_likelihoods[ref_idx];
        let la = row.haplotype_log10_likelihoods[alt_idx];
        if lr >= la {
            ref_ad += 1;
        } else {
            alt_ad += 1;
        }
    }
    if ref_ad == 0 && alt_ad == 0 {
        ref_ad = 1;
    }
    let gls = if legacy_gl {
        biallelic_genotype_log10_likelihoods_parity_legacy(rows, ref_idx, alt_idx)
    } else {
        biallelic_genotype_log10_likelihoods_gatk(rows, ref_idx, alt_idx)
    };
    let format = emit_genotype_format_fields(&gls, &[ref_ad, alt_ad])?;
    Ok(ParityRegionGenotypeDump {
        haplotype_count: hap_count,
        read_count,
        genotyped: true,
        ref_hap_idx: ref_idx,
        alt_hap_idx: alt_idx,
        best_hap_idx: best,
        genotype_log10: gls,
        format: Some(format),
    })
}

pub fn write_parity_region_genotype_dump(
    out: &mut impl Write,
    dump: &ParityRegionGenotypeDump,
) -> GatkResult<()> {
    writeln!(out, "haplotype_count\t{}", dump.haplotype_count)
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    writeln!(out, "read_count\t{}", dump.read_count)
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    writeln!(out, "genotyped\t{}", dump.genotyped)
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    let Some(format) = &dump.format else {
        return Ok(());
    };
    writeln!(out, "ref_haplotype_index\t{}", dump.ref_hap_idx)
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    writeln!(out, "alt_haplotype_index\t{}", dump.alt_hap_idx)
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    writeln!(out, "best_haplotype_index\t{}", dump.best_hap_idx)
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    for (i, gl) in dump.genotype_log10.iter().enumerate() {
        writeln!(out, "genotype_{i}_log10\t{gl}")
            .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    }
    writeln!(
        out,
        "pl\t{}",
        format
            .pl
            .iter()
            .map(|p| p.to_string())
            .collect::<Vec<_>>()
            .join(",")
    )
    .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    writeln!(out, "gq\t{}", format.gq)
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    writeln!(
        out,
        "ad\t{}",
        format
            .ad
            .iter()
            .map(|a| a.to_string())
            .collect::<Vec<_>>()
            .join(",")
    )
    .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    writeln!(out, "dp\t{}", format.dp)
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn parity_region_genotype_dump_writes_genotyped_row() {
        let rows = vec![ReadLikelihoodRow {
            read_id: "r1".into(),
            haplotype_log10_likelihoods: vec![-1.0],
        }];
        let dump = parity_region_genotype_from_rows(&rows, &[true]).expect("dump");
        assert!(dump.genotyped);
        let mut buf = Vec::new();
        write_parity_region_genotype_dump(&mut buf, &dump).expect("write");
        let s = String::from_utf8(buf).expect("utf8");
        assert!(s.contains("genotyped\ttrue"), "output:\n{s}");
    }
}
