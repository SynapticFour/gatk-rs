//! GenotypeGVCFs — joint genotyping of a multi-sample gVCF (post CombineGVCFs).
//! Observable Java contract (GATK 4.4 `GenotypeGVCFs` / `GenotypeGVCFsEngine`):
//! Input: one multi-sample gVCF with diploid PL/AD (typically CombineGVCFs output).
//! Drop symbolic `<NON_REF>` before genotyping; remap PL/AD onto remaining alleles.
//! Estimate cohort allele frequencies via Dirichlet / HWE EM over all sample GLs
//! ([`crate::af_calc`]), derive site QUAL = −10·log10 P(no variant), assign per-sample
//! genotypes from PL + AF prior (preferring PLs), emit AC/AN/AF/NS/DP/QD.
//! Default emit threshold: `standard_confidence_for_calling` = 30.0 (hom-ref blocks omitted).
//! Annotation reuse: core INFO via [`crate::genotyping::compute_core_variant_annotations`];
//! QD = QUAL / variant depth. Read-based FS/MQ require pileups and are not computed here.

use crate::af_calc::{
    calculate_multiallelic_af_em, qual_from_log10_p_no_variant, AfCalculatorConfig,
};
use crate::genotyping::{
    compute_core_variant_annotations, gq_phred_from_pl, SampleAnnotationInput,
};
use crate::ref_confidence_merger::NON_REF_ALLELE;
use gatk_common::{GatkError, GatkResult};
use gatk_core::io::vcf::{
    Contig, FormatField, Genotype, InfoField, InfoValue, SampleData, VcfHeader, VcfReader,
    VcfRecord, VcfWriter,
};
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use std::collections::HashMap;
use std::path::{Path, PathBuf};

/// Default GATK 4.1+ GenotypeGVCFs call confidence.
pub const DEFAULT_STAND_CALL_CONF: f64 = 30.0;

/// CLI / library options for GenotypeGVCFs.
#[derive(Debug, Clone)]
pub struct GenotypeGvcfsArgs {
    pub reference: PathBuf,
    pub variant: PathBuf,
    pub output: PathBuf,
    pub intervals: Option<String>,
    /// Minimum QUAL to emit a site (GATK `-stand-call-conf`).
    pub stand_call_conf: f64,
    /// When true, emit sites that fail the call confidence (for debugging).
    pub include_non_variant_sites: bool,
}

impl Default for GenotypeGvcfsArgs {
    fn default() -> Self {
        Self {
            reference: PathBuf::new(),
            variant: PathBuf::new(),
            output: PathBuf::new(),
            intervals: None,
            stand_call_conf: DEFAULT_STAND_CALL_CONF,
            include_non_variant_sites: false,
        }
    }
}

/// Result of genotyping one gVCF site (tests / harness).
#[derive(Debug, Clone)]
pub struct GenotypedSite {
    pub record: VcfRecord,
    pub qual: f64,
    pub allele_frequencies: Vec<f64>,
    pub log10_p_no_variant: f64,
}

/// Run GenotypeGVCFs from filesystem paths (CLI entry).
pub fn run_genotype_gvcfs(args: &GenotypeGvcfsArgs) -> GatkResult<()> {
    let dictionary = SequenceDictionary::from_fasta_path(&args.reference)?;
    let clip = resolve_clip(&dictionary, args.intervals.as_deref())?;
    let mut reader = VcfReader::from_file(&args.variant)?;
    let in_header = reader.header().clone();
    let samples = in_header.samples.clone();
    let raw = reader.read_all_records()?;
    let genotyped = genotype_gvcf_records(
        &raw,
        &samples,
        &GenotypeGvcfsConfig {
            stand_call_conf: args.stand_call_conf,
            include_non_variant_sites: args.include_non_variant_sites,
            clip,
        },
    )?;
    let out_recs: Vec<VcfRecord> = genotyped.into_iter().map(|g| g.record).collect();
    let header = build_output_header(&in_header, Some(args.reference.as_path()))?;
    let mut writer = VcfWriter::new(&args.output, header)?;
    writer.write_header()?;
    writer.write_records(&out_recs)?;
    Ok(())
}

/// Tunables for in-memory genotyping.
#[derive(Debug, Clone, Default)]
pub struct GenotypeGvcfsConfig {
    pub stand_call_conf: f64,
    pub include_non_variant_sites: bool,
    pub clip: Option<HashMap<String, (u64, u64)>>,
}

impl GenotypeGvcfsConfig {
    pub fn with_defaults() -> Self {
        Self {
            stand_call_conf: DEFAULT_STAND_CALL_CONF,
            include_non_variant_sites: false,
            clip: None,
        }
    }
}

/// Joint-genotype a slice of multi-sample gVCF records.
pub fn genotype_gvcf_records(
    records: &[VcfRecord],
    sample_names: &[String],
    config: &GenotypeGvcfsConfig,
) -> GatkResult<Vec<GenotypedSite>> {
    let mut out = Vec::new();
    for rec in records {
        if let Some(clip) = &config.clip {
            if let Some(&(lo, hi)) = clip.get(&rec.chromosome) {
                if rec.position < lo || rec.position > hi {
                    continue;
                }
            } else {
                continue;
            }
        }
        if let Some(site) = genotype_one_record(rec, sample_names, config)? {
            out.push(site);
        }
    }
    Ok(out)
}

fn resolve_clip(
    dictionary: &SequenceDictionary,
    intervals: Option<&str>,
) -> GatkResult<Option<HashMap<String, (u64, u64)>>> {
    let Some(iv) = intervals else {
        return Ok(None);
    };
    let specs = parse_intervals_cli_string(dictionary, iv)?;
    let mut map = HashMap::new();
    for spec in specs {
        let start = spec.start.unwrap_or(1);
        let end = match spec.end {
            Some(e) => e,
            None => dictionary
                .contig(&spec.contig)
                .map(|c| c.length)
                .ok_or_else(|| GatkError::argument(format!("unknown contig {}", spec.contig)))?,
        };
        map.insert(spec.contig, (start, end));
    }
    Ok(Some(map))
}

fn genotype_one_record(
    rec: &VcfRecord,
    sample_names: &[String],
    config: &GenotypeGvcfsConfig,
) -> GatkResult<Option<GenotypedSite>> {
    let (keep_alleles, allele_map) = alleles_without_non_ref(&rec.reference, &rec.alternate);
    // keep_alleles[0] = REF; rest = concrete ALTs (may be empty → hom-ref block).
    if keep_alleles.len() < 2 {
        // Pure reference-confidence block — GenotypeGVCFs omits by default.
        if !config.include_non_variant_sites {
            return Ok(None);
        }
    }

    let n_alleles = keep_alleles.len();
    let n_alt = n_alleles.saturating_sub(1);
    if n_alt == 0 {
        return Ok(None);
    }

    let old_n = 1 + rec.alternate.len();
    let mut sample_gls: Vec<Vec<f64>> = Vec::with_capacity(rec.samples.len());
    let mut remapped_pls: Vec<Option<Vec<i32>>> = Vec::with_capacity(rec.samples.len());
    let mut remapped_ads: Vec<Option<Vec<u32>>> = Vec::with_capacity(rec.samples.len());

    for sample in &rec.samples {
        let pl_i32 = sample.pl.as_ref().map(|pl| {
            pl.iter()
                .map(|&x| i32::try_from(x).unwrap_or(i32::MAX))
                .collect::<Vec<_>>()
        });
        let new_pl = pl_i32
            .as_ref()
            .and_then(|pl| remap_pl_drop_alleles(pl, old_n, &allele_map, n_alleles));
        let new_ad = sample
            .ad
            .as_ref()
            .map(|ad| remap_ad(ad, &allele_map, n_alleles));
        let gl = new_pl
            .as_ref()
            .map(|pl| pl.iter().map(|&p| -(p as f64) / 10.0).collect())
            .unwrap_or_default();
        remapped_pls.push(new_pl);
        remapped_ads.push(new_ad);
        sample_gls.push(gl);
    }

    let gl_refs: Vec<&[f64]> = sample_gls
        .iter()
        .filter(|g| g.len() == n_alleles * (n_alleles + 1) / 2)
        .map(|g| g.as_slice())
        .collect();

    if gl_refs.is_empty() {
        return Ok(None);
    }

    let af = calculate_multiallelic_af_em(&gl_refs, n_alleles, &AfCalculatorConfig::default())?;
    let qual = qual_from_log10_p_no_variant(af.log10_posterior_no_variant);

    if qual < config.stand_call_conf && !config.include_non_variant_sites {
        return Ok(None);
    }

    let pairs = diploid_pairs(n_alleles);
    let log10_af: Vec<f64> = af
        .allele_frequencies
        .iter()
        .map(|&f| f.max(1e-300).log10())
        .collect();

    let mut out_samples = Vec::with_capacity(rec.samples.len());
    let mut ann_inputs = Vec::with_capacity(rec.samples.len());
    let mut qd_depth = 0_i32;

    for (idx, sample) in rec.samples.iter().enumerate() {
        // Take ownership — each remapped vector is consumed once into the output sample.
        let pl = std::mem::take(&mut remapped_pls[idx]);
        let ad = std::mem::take(&mut remapped_ads[idx]);
        let (gt, gq) = if let Some(ref plv) = pl {
            if sample_gls[idx].len() == pairs.len() {
                let post = genotype_posteriors_hwe(&sample_gls[idx], &log10_af, &pairs);
                let best = post
                    .iter()
                    .enumerate()
                    .max_by(|a, b| a.1.total_cmp(b.1))
                    .map(|(i, _)| i)
                    .unwrap_or(0);
                // Prefer PL min when it agrees with MAP; else use MAP under AF prior.
                let pl_best = plv
                    .iter()
                    .enumerate()
                    .min_by_key(|(_, &p)| p)
                    .map(|(i, _)| i)
                    .unwrap_or(best);
                let chosen = if pl_best == best { pl_best } else { best };
                let (a, b) = pairs[chosen];
                let gq = gq_phred_from_pl(plv);
                (
                    Some(Genotype {
                        alleles: vec![a as i32, b as i32],
                        phased: false,
                    }),
                    Some(gq as f64),
                )
            } else {
                (None, None)
            }
        } else {
            (None, None)
        };

        if let Some(ref g) = gt {
            if g.alleles.iter().any(|&a| a > 0) {
                if let Some(ref adv) = ad {
                    let total: i32 = adv.iter().map(|&x| x as i32).sum();
                    qd_depth += total.max(0);
                } else if let Some(dp) = sample.dp {
                    qd_depth += dp as i32;
                }
            }
            ann_inputs.push(SampleAnnotationInput {
                // CLONE: needed — annotation input owns diploid allele indices (cheap `Vec<i32>` of len 2).
                genotype_alleles: g.alleles.clone(),
                dp: sample.dp.map(|d| d as i32),
            });
        } else {
            ann_inputs.push(SampleAnnotationInput {
                genotype_alleles: vec![-1, -1],
                dp: sample.dp.map(|d| d as i32),
            });
        }

        let pl_u32 = pl.map(|v| {
            v.into_iter()
                .map(|x| u32::try_from(x.max(0)).unwrap_or(u32::MAX))
                .collect()
        });
        out_samples.push(SampleData {
            gt,
            gq,
            dp: sample.dp,
            ad,
            pl: pl_u32,
            other: Vec::new(),
        });
    }

    let core = compute_core_variant_annotations(n_alt, &ann_inputs)?;
    let qd = if qd_depth > 0 {
        Some(qual / qd_depth as f64)
    } else {
        None
    };

    let alts: Vec<String> = keep_alleles.iter().skip(1).cloned().collect();
    let mut info = vec![
        InfoValue::Integer("AC".to_string(), core.ac.clone()),
        InfoValue::Integer("AN".to_string(), vec![core.an]),
        InfoValue::Float("AF".to_string(), core.af.clone()),
        InfoValue::Integer("NS".to_string(), vec![core.ns]),
        InfoValue::Integer("DP".to_string(), vec![core.dp]),
    ];
    if let Some(qd_v) = qd {
        info.push(InfoValue::Float("QD".to_string(), vec![qd_v]));
    }
    // MLEAC / MLEAF from AF calculator expected counts.
    if af.mle_allele_counts.len() == n_alleles {
        let mleac: Vec<i32> = af.mle_allele_counts[1..]
            .iter()
            .map(|c| c.round() as i32)
            .collect();
        let total: f64 = af.mle_allele_counts.iter().sum();
        let mleaf: Vec<f64> = if total > 0.0 {
            af.mle_allele_counts[1..]
                .iter()
                .map(|c| c / total)
                .collect()
        } else {
            vec![0.0; n_alt]
        };
        info.push(InfoValue::Integer("MLEAC".to_string(), mleac));
        info.push(InfoValue::Float("MLEAF".to_string(), mleaf));
    }

    if !sample_names.is_empty() && sample_names.len() != out_samples.len() {
        return Err(GatkError::argument(format!(
            "GenotypeGVCFs sample count mismatch: header {} vs record {}",
            sample_names.len(),
            out_samples.len()
        )));
    }
    Ok(Some(GenotypedSite {
        record: VcfRecord {
            chromosome: rec.chromosome.clone(),
            position: rec.position,
            id: rec.id.clone(),
            reference: keep_alleles[0].clone(),
            alternate: alts,
            quality: Some(qual),
            filter: vec![".".to_string()],
            info,
            format: vec![
                "GT".to_string(),
                "GQ".to_string(),
                "DP".to_string(),
                "AD".to_string(),
                "PL".to_string(),
            ],
            samples: out_samples,
        },
        qual,
        allele_frequencies: af.allele_frequencies,
        log10_p_no_variant: af.log10_posterior_no_variant,
    }))
}

fn alleles_without_non_ref(reference: &str, alts: &[String]) -> (Vec<String>, Vec<Option<usize>>) {
    let mut keep = vec![reference.to_string()];
    let mut map = vec![Some(0usize)]; // old REF → new 0
    for alt in alts {
        if alt == NON_REF_ALLELE {
            map.push(None);
        } else {
            map.push(Some(keep.len()));
            // CLONE: needed — `alts` is borrowed; keep-list owns concrete allele strings.
            keep.push(alt.clone());
        }
    }
    (keep, map)
}

fn diploid_pairs(n_alleles: usize) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();
    for j in 0..n_alleles {
        for i in 0..=j {
            pairs.push((i, j));
        }
    }
    pairs
}

fn remap_pl_drop_alleles(
    old_pl: &[i32],
    old_n: usize,
    allele_map: &[Option<usize>],
    new_n: usize,
) -> Option<Vec<i32>> {
    let old_pairs = diploid_pairs(old_n);
    let new_pairs = diploid_pairs(new_n);
    if old_pl.len() != old_pairs.len() {
        return None;
    }
    // For each new genotype, find best (min PL) among old genotypes whose alleles
    // map into the new pair (NON_REF contributions collapse to nearest kept allele).
    let mut out = Vec::with_capacity(new_pairs.len());
    for &(ni, nj) in &new_pairs {
        let mut best = i32::MAX;
        for (oi, &(oa, ob)) in old_pairs.iter().enumerate() {
            let ma = allele_map.get(oa).copied().flatten();
            let mb = allele_map.get(ob).copied().flatten();
            match (ma, mb) {
                (Some(a), Some(b)) => {
                    let (a, b) = if a <= b { (a, b) } else { (b, a) };
                    if a == ni && b == nj {
                        best = best.min(old_pl[oi]);
                    }
                }
                // Genotypes involving dropped NON_REF: fold into hom-ref / het involving mapped allele.
                (Some(a), None) | (None, Some(a)) => {
                    if ni == a && nj == a {
                        best = best.min(old_pl[oi]);
                    } else if (ni == 0 && nj == a) || (ni == a && nj == 0) {
                        best = best.min(old_pl[oi]);
                    }
                }
                (None, None) => {
                    if ni == 0 && nj == 0 {
                        best = best.min(old_pl[oi]);
                    }
                }
            }
        }
        if best == i32::MAX {
            // Fallback: use NON_REF-mapped index like AlleleSubsetting (to NON_REF→ absent → max PL).
            best = old_pl.iter().copied().max().unwrap_or(0);
        }
        out.push(best);
    }
    // Renormalize so min PL = 0.
    if let Some(&m) = out.iter().min() {
        for p in &mut out {
            *p -= m;
        }
    }
    Some(out)
}

fn remap_ad(old_ad: &[u32], allele_map: &[Option<usize>], new_n: usize) -> Vec<u32> {
    let mut out = vec![0u32; new_n];
    for (old_i, &count) in old_ad.iter().enumerate() {
        if let Some(Some(new_i)) = allele_map.get(old_i).copied() {
            if new_i < new_n {
                out[new_i] = out[new_i].saturating_add(count);
            }
        }
    }
    out
}

fn genotype_posteriors_hwe(
    log10_likelihoods: &[f64],
    log10_af: &[f64],
    pairs: &[(usize, usize)],
) -> Vec<f64> {
    use crate::activity_scoring::log10_sum_log10;
    let n = pairs.len().min(log10_likelihoods.len());
    let mut log10_post = vec![f64::NEG_INFINITY; n];
    for (gi, &(a, b)) in pairs.iter().enumerate().take(n) {
        let log10_hwe = if a == b {
            2.0 * log10_af[a]
        } else {
            (2.0_f64).log10() + log10_af[a] + log10_af[b]
        };
        log10_post[gi] = log10_likelihoods[gi] + log10_hwe;
    }
    let s = log10_sum_log10(&log10_post);
    log10_post.into_iter().map(|x| 10_f64.powf(x - s)).collect()
}

fn build_output_header(input: &VcfHeader, reference: Option<&Path>) -> GatkResult<VcfHeader> {
    let mut header = VcfHeader::default();
    header.file_format = "VCFv4.2".to_string();
    header.source = Some("gatk-rs GenotypeGVCFs".to_string());
    if let Some(r) = reference {
        header.reference = Some(r.display().to_string());
    }
    header.contigs = input
        .contigs
        .iter()
        .filter(|c| !c.id.is_empty())
        .cloned()
        .collect();
    if header.contigs.is_empty() {
        header.contigs.push(Contig {
            id: "chr1".to_string(),
            length: None,
            md5: None,
            assembly: None,
            species: None,
            uri: None,
        });
    }
    header.info_fields = vec![
        info("AC", "A", "Integer", "Allele count in genotypes"),
        info("AF", "A", "Float", "Allele frequency"),
        info(
            "AN",
            "1",
            "Integer",
            "Total number of alleles in called genotypes",
        ),
        info("DP", "1", "Integer", "Approximate read depth"),
        info("NS", "1", "Integer", "Number of samples with data"),
        info("QD", "1", "Float", "Variant Confidence/Quality by Depth"),
        info(
            "MLEAC",
            "A",
            "Integer",
            "Maximum likelihood expectation (MLE) for the allele counts",
        ),
        info(
            "MLEAF",
            "A",
            "Float",
            "Maximum likelihood expectation (MLE) for the allele frequency",
        ),
    ];
    header.format_fields = vec![
        fmt("GT", "1", "String", "Genotype"),
        fmt("GQ", "1", "Integer", "Genotype Quality"),
        fmt("DP", "1", "Integer", "Approximate read depth"),
        fmt("AD", "R", "Integer", "Allelic depths"),
        fmt(
            "PL",
            "G",
            "Integer",
            "Normalized, Phred-scaled likelihoods for genotypes",
        ),
    ];
    header.samples = input.samples.clone();
    Ok(header)
}

fn info(id: &str, number: &str, ty: &str, desc: &str) -> InfoField {
    InfoField {
        id: id.to_string(),
        number: number.to_string(),
        type_field: ty.to_string(),
        description: desc.to_string(),
        source: None,
        version: None,
    }
}

fn fmt(id: &str, number: &str, ty: &str, desc: &str) -> FormatField {
    FormatField {
        id: id.to_string(),
        number: number.to_string(),
        type_field: ty.to_string(),
        description: desc.to_string(),
    }
}

// --------------------------------------------------------------------------
// Tests
// --------------------------------------------------------------------------

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;

    fn site(pos: u64, refr: &str, alts: &[&str], samples: Vec<SampleData>) -> VcfRecord {
        VcfRecord {
            chromosome: "chr1".to_string(),
            position: pos,
            id: ".".to_string(),
            reference: refr.to_string(),
            alternate: alts.iter().map(|s| s.to_string()).collect(),
            quality: None,
            filter: vec![".".to_string()],
            info: Vec::new(),
            format: vec![
                "GT".to_string(),
                "PL".to_string(),
                "DP".to_string(),
                "AD".to_string(),
            ],
            samples,
        }
    }

    fn sample_pl(pl: &[u32], dp: u32, ad: &[u32]) -> SampleData {
        SampleData {
            gt: None,
            gq: None,
            dp: Some(dp),
            ad: Some(ad.to_vec()),
            pl: Some(pl.to_vec()),
            other: Vec::new(),
        }
    }

    #[test]
    fn t01_biallelic_het_cohort_calls_variant() {
        // 3 samples: two hets, one hom-ref — expect emit with AF≈1/3.
        let names = vec!["A".into(), "B".into(), "C".into()];
        let rec = site(
            100,
            "A",
            &["G", NON_REF_ALLELE],
            vec![
                sample_pl(&[100, 0, 100, 100, 100, 100], 20, &[10, 10, 0]),
                sample_pl(&[100, 0, 100, 100, 100, 100], 20, &[10, 10, 0]),
                sample_pl(&[0, 100, 1000, 100, 1000, 1000], 20, &[20, 0, 0]),
            ],
        );
        let mut cfg = GenotypeGvcfsConfig::with_defaults();
        cfg.stand_call_conf = 0.0; // always emit for test
        let out = genotype_gvcf_records(&[rec], &names, &cfg).unwrap();
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].record.alternate, vec!["G".to_string()]);
        assert!(out[0].qual > 0.0);
        assert!(out[0].allele_frequencies[1] > 0.1);
        assert!(out[0].allele_frequencies[1] < 0.6);
    }

    #[test]
    fn t02_homref_block_skipped_by_default() {
        let names = vec!["A".into(), "B".into()];
        let rec = site(
            1,
            "A",
            &[NON_REF_ALLELE],
            vec![
                sample_pl(&[0, 90, 900], 10, &[10, 0]),
                sample_pl(&[0, 60, 600], 10, &[10, 0]),
            ],
        );
        let out =
            genotype_gvcf_records(&[rec], &names, &GenotypeGvcfsConfig::with_defaults()).unwrap();
        assert!(out.is_empty());
    }

    #[test]
    fn t03_stand_call_conf_filters_weak_sites() {
        let names = vec!["A".into()];
        // Ambiguous PL → low QUAL.
        let rec = site(
            50,
            "C",
            &["T", NON_REF_ALLELE],
            vec![sample_pl(&[10, 0, 10, 10, 10, 10], 4, &[2, 2, 0])],
        );
        let mut cfg = GenotypeGvcfsConfig::with_defaults();
        cfg.stand_call_conf = 1000.0;
        let out = genotype_gvcf_records(&[rec], &names, &cfg).unwrap();
        assert!(out.is_empty());
    }

    #[test]
    fn t04_multiallelic_keeps_both_alts() {
        let names = vec!["S1".into(), "S2".into()];
        let rec = site(
            10,
            "A",
            &["G", "T", NON_REF_ALLELE],
            vec![
                // 10 PLs for 4 alleles
                sample_pl(
                    &[100, 0, 100, 100, 100, 100, 100, 100, 100, 100],
                    20,
                    &[10, 10, 0, 0],
                ),
                sample_pl(&[90, 90, 90, 0, 90, 90, 90, 90, 90, 90], 16, &[8, 0, 8, 0]),
            ],
        );
        let mut cfg = GenotypeGvcfsConfig::with_defaults();
        cfg.stand_call_conf = 0.0;
        let out = genotype_gvcf_records(&[rec], &names, &cfg).unwrap();
        assert_eq!(out.len(), 1);
        assert!(out[0].record.alternate.contains(&"G".to_string()));
        assert!(out[0].record.alternate.contains(&"T".to_string()));
        assert!(!out[0].record.alternate.iter().any(|a| a == NON_REF_ALLELE));
    }

    #[test]
    fn t05_non_ref_stripped_from_pl_length() {
        let names = vec!["A".into()];
        let rec = site(
            1,
            "A",
            &["C", NON_REF_ALLELE],
            vec![sample_pl(&[0, 100, 1000, 100, 1000, 1000], 10, &[10, 0, 0])],
        );
        let mut cfg = GenotypeGvcfsConfig::with_defaults();
        cfg.stand_call_conf = 0.0;
        let out = genotype_gvcf_records(&[rec], &names, &cfg).unwrap();
        // After drop NON_REF: 3 PLs
        assert_eq!(out[0].record.samples[0].pl.as_ref().unwrap().len(), 3);
    }

    #[test]
    fn t06_ac_an_af_annotations() {
        let names = vec!["A".into(), "B".into(), "C".into()];
        let rec = site(
            7,
            "G",
            &["A", NON_REF_ALLELE],
            vec![
                sample_pl(&[1000, 0, 1000, 1000, 1000, 1000], 30, &[15, 15, 0]),
                sample_pl(&[1000, 0, 1000, 1000, 1000, 1000], 30, &[15, 15, 0]),
                sample_pl(&[1000, 0, 1000, 1000, 1000, 1000], 30, &[15, 15, 0]),
            ],
        );
        let mut cfg = GenotypeGvcfsConfig::with_defaults();
        cfg.stand_call_conf = 0.0;
        let out = genotype_gvcf_records(&[rec], &names, &cfg).unwrap();
        let ac = out[0]
            .record
            .info
            .iter()
            .find_map(|i| match i {
                InfoValue::Integer(id, v) if id == "AC" => Some(v[0]),
                _ => None,
            })
            .unwrap();
        let an = out[0]
            .record
            .info
            .iter()
            .find_map(|i| match i {
                InfoValue::Integer(id, v) if id == "AN" => Some(v[0]),
                _ => None,
            })
            .unwrap();
        assert_eq!(an, 6);
        assert_eq!(ac, 3); // three hets → 3 alt alleles
    }

    #[test]
    fn t07_missing_pl_sample_is_nocall() {
        let names = vec!["A".into(), "B".into()];
        let rec = site(
            3,
            "T",
            &["C", NON_REF_ALLELE],
            vec![
                sample_pl(&[100, 0, 100, 100, 100, 100], 20, &[10, 10, 0]),
                SampleData {
                    gt: None,
                    gq: None,
                    dp: None,
                    ad: None,
                    pl: None,
                    other: Vec::new(),
                },
            ],
        );
        let mut cfg = GenotypeGvcfsConfig::with_defaults();
        cfg.stand_call_conf = 0.0;
        let out = genotype_gvcf_records(&[rec], &names, &cfg).unwrap();
        assert!(out[0].record.samples[1].gt.is_none());
        assert!(out[0].record.samples[0].gt.is_some());
    }

    #[test]
    fn t08_ten_sample_cohort_af_near_half() {
        // 5 het + 5 hom-ref → AF ≈ 0.25
        let mut names = Vec::new();
        let mut samples = Vec::new();
        for i in 0..5 {
            names.push(format!("H{i}"));
            samples.push(sample_pl(&[200, 0, 200, 200, 200, 200], 20, &[10, 10, 0]));
        }
        for i in 0..5 {
            names.push(format!("R{i}"));
            samples.push(sample_pl(&[0, 200, 2000, 200, 2000, 2000], 20, &[20, 0, 0]));
        }
        let rec = site(99, "A", &["G", NON_REF_ALLELE], samples);
        let mut cfg = GenotypeGvcfsConfig::with_defaults();
        cfg.stand_call_conf = 0.0;
        let out = genotype_gvcf_records(&[rec], &names, &cfg).unwrap();
        let af = out[0].allele_frequencies[1];
        assert!(af > 0.15 && af < 0.40, "af={af}");
    }

    #[test]
    fn t09_qd_present_when_variant_depth() {
        let names = vec!["A".into()];
        let rec = site(
            1,
            "A",
            &["T", NON_REF_ALLELE],
            vec![sample_pl(&[500, 0, 500, 500, 500, 500], 40, &[20, 20, 0])],
        );
        let mut cfg = GenotypeGvcfsConfig::with_defaults();
        cfg.stand_call_conf = 0.0;
        let out = genotype_gvcf_records(&[rec], &names, &cfg).unwrap();
        assert!(out[0]
            .record
            .info
            .iter()
            .any(|i| matches!(i, InfoValue::Float(id, _) if id == "QD")));
    }

    #[test]
    fn t10_hom_var_cohort_high_af() {
        let names = vec!["A".into(), "B".into(), "C".into()];
        let rec = site(
            2,
            "C",
            &["G", NON_REF_ALLELE],
            vec![
                sample_pl(&[1000, 1000, 0, 1000, 1000, 1000], 30, &[0, 30, 0]),
                sample_pl(&[1000, 1000, 0, 1000, 1000, 1000], 30, &[0, 30, 0]),
                sample_pl(&[1000, 1000, 0, 1000, 1000, 1000], 30, &[0, 30, 0]),
            ],
        );
        let mut cfg = GenotypeGvcfsConfig::with_defaults();
        cfg.stand_call_conf = 0.0;
        let out = genotype_gvcf_records(&[rec], &names, &cfg).unwrap();
        assert!(out[0].allele_frequencies[1] > 0.8);
        for s in &out[0].record.samples {
            assert_eq!(s.gt.as_ref().unwrap().alleles, vec![1, 1]);
        }
    }

    #[test]
    fn t11_clip_interval_filters_positions() {
        let names = vec!["A".into()];
        let r1 = site(
            10,
            "A",
            &["G", NON_REF_ALLELE],
            vec![sample_pl(&[100, 0, 100, 100, 100, 100], 10, &[5, 5, 0])],
        );
        let r2 = site(
            50,
            "A",
            &["G", NON_REF_ALLELE],
            vec![sample_pl(&[100, 0, 100, 100, 100, 100], 10, &[5, 5, 0])],
        );
        let mut cfg = GenotypeGvcfsConfig::with_defaults();
        cfg.stand_call_conf = 0.0;
        let mut clip = HashMap::new();
        clip.insert("chr1".to_string(), (40, 60));
        cfg.clip = Some(clip);
        let out = genotype_gvcf_records(&[r1, r2], &names, &cfg).unwrap();
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].record.position, 50);
    }

    #[test]
    fn t12_qual_increases_with_more_carriers() {
        let mk = |n_het: usize| {
            let mut names = Vec::new();
            let mut samples = Vec::new();
            for i in 0..n_het {
                names.push(format!("H{i}"));
                samples.push(sample_pl(&[300, 0, 300, 300, 300, 300], 20, &[10, 10, 0]));
            }
            for i in 0..(8 - n_het) {
                names.push(format!("R{i}"));
                samples.push(sample_pl(&[0, 300, 3000, 300, 3000, 3000], 20, &[20, 0, 0]));
            }
            let rec = site(1, "A", &["C", NON_REF_ALLELE], samples);
            let mut cfg = GenotypeGvcfsConfig::with_defaults();
            cfg.stand_call_conf = 0.0;
            genotype_gvcf_records(&[rec], &names, &cfg).unwrap()[0].qual
        };
        assert!(mk(4) > mk(1));
    }

    #[test]
    fn t13_gq_capped_at_99() {
        let names = vec!["A".into()];
        let rec = site(
            1,
            "A",
            &["G", NON_REF_ALLELE],
            vec![sample_pl(
                &[5000, 0, 5000, 5000, 5000, 5000],
                50,
                &[25, 25, 0],
            )],
        );
        let mut cfg = GenotypeGvcfsConfig::with_defaults();
        cfg.stand_call_conf = 0.0;
        let out = genotype_gvcf_records(&[rec], &names, &cfg).unwrap();
        assert!(out[0].record.samples[0].gq.unwrap() <= 99.0);
    }

    #[test]
    fn t14_empty_input_ok() {
        let out = genotype_gvcf_records(&[], &[], &GenotypeGvcfsConfig::with_defaults()).unwrap();
        assert!(out.is_empty());
    }

    #[test]
    fn t15_mleaf_sums_near_alt_af() {
        let names = vec!["A".into(), "B".into()];
        let rec = site(
            1,
            "A",
            &["G", NON_REF_ALLELE],
            vec![
                sample_pl(&[100, 0, 100, 100, 100, 100], 20, &[10, 10, 0]),
                sample_pl(&[0, 100, 1000, 100, 1000, 1000], 20, &[20, 0, 0]),
            ],
        );
        let mut cfg = GenotypeGvcfsConfig::with_defaults();
        cfg.stand_call_conf = 0.0;
        let out = genotype_gvcf_records(&[rec], &names, &cfg).unwrap();
        let mleaf = out[0]
            .record
            .info
            .iter()
            .find_map(|i| match i {
                InfoValue::Float(id, v) if id == "MLEAF" => Some(v[0]),
                _ => None,
            })
            .unwrap();
        assert!((mleaf - out[0].allele_frequencies[1]).abs() < 0.25);
    }

    #[test]
    fn t16_allele_map_drops_only_non_ref() {
        let (keep, map) =
            alleles_without_non_ref("A", &["G".into(), NON_REF_ALLELE.into(), "T".into()]);
        assert_eq!(keep, vec!["A", "G", "T"]);
        assert_eq!(map, vec![Some(0), Some(1), None, Some(2)]);
    }
}
