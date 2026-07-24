//! CombineGVCFs — merge per-sample gVCFs into a multi-sample gVCF (no joint genotyping).
//! Observable Java contract (GATK 4.4 `CombineGVCFs`):
//! Union sample columns from input gVCFs.
//! At each genomic breakpoint induced by record starts / `END`+1, merge overlapping
//! reference-confidence VCs via `ReferenceConfidenceVariantContextMerger` semantics
//! ([`crate::ref_confidence_merger::merge_reference_confidence`]).
//! Remap PL/AD onto the unified allele list; missing samples emit no-call genotypes.
//! Hom-ref `<NON_REF>` blocks carry `INFO/END`; variant sites emit discrete rows.
//! This module is a Rust-native walker over that merger kernel — not a Java class clone.

use crate::ref_confidence_merger::{
    merge_reference_confidence, MergeAllele, MergeGenotype, MergeVcInput, NON_REF_ALLELE,
};
use gatk_common::{GatkError, GatkResult};
use gatk_core::io::vcf::{
    Contig, FormatField, Genotype, InfoField, InfoValue, SampleData, VcfHeader, VcfReader,
    VcfRecord, VcfWriter,
};
use gatk_core::reference::{reference_base_at_1based, SequenceDictionary};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};

/// CLI / library options for CombineGVCFs.
#[derive(Debug, Clone)]
pub struct CombineGvcfsArgs {
    pub reference: PathBuf,
    pub variant_paths: Vec<PathBuf>,
    pub output: PathBuf,
    /// Optional `-L` interval string (same grammar as HaplotypeCaller).
    pub intervals: Option<String>,
}

/// One loaded single-sample gVCF site (variant or reference block).
#[derive(Debug, Clone)]
pub struct GvcfSite {
    pub contig: String,
    /// 1-based inclusive start (VCF POS).
    pub start: u64,
    /// 1-based inclusive end (`INFO/END` or POS + REF.len - 1).
    pub end: u64,
    pub ref_allele: String,
    pub alts: Vec<String>,
    pub sample: SampleData,
    pub source: String,
}

/// One input gVCF after load.
#[derive(Debug, Clone)]
pub struct LoadedGvcf {
    pub path: PathBuf,
    pub sample_name: String,
    pub header: VcfHeader,
    pub sites: Vec<GvcfSite>,
}

/// Inclusive genomic interval on one contig.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Span {
    start: u64,
    end: u64,
}

/// Run CombineGVCFs from filesystem paths (CLI entry).
pub fn run_combine_gvcfs(args: &CombineGvcfsArgs) -> GatkResult<()> {
    if args.variant_paths.is_empty() {
        return Err(GatkError::argument(
            "CombineGVCFs requires at least one -V input",
        ));
    }
    let dictionary = SequenceDictionary::from_fasta_path(&args.reference)?;
    let clip = resolve_clip_spans(&dictionary, args.intervals.as_deref())?;
    let mut loaded = Vec::with_capacity(args.variant_paths.len());
    for path in &args.variant_paths {
        loaded.push(load_gvcf(path)?);
    }
    let records = combine_loaded_gvcfs(&loaded, Some(args.reference.as_path()), clip.as_ref())?;
    let header = build_combined_header(&loaded, Some(args.reference.as_path()))?;
    let mut writer = VcfWriter::new(&args.output, header)?;
    writer.write_header()?;
    writer.write_records(&records)?;
    Ok(())
}

/// Merge already-loaded gVCFs (unit-test / harness entry).
/// `reference` supplies bases when no input VC starts at a locus (spanning-only merges).
/// `clip` restricts emission to inclusive spans per contig when present.
pub fn combine_loaded_gvcfs(
    inputs: &[LoadedGvcf],
    reference: Option<&Path>,
    clip: Option<&HashMap<String, Vec<Span>>>,
) -> GatkResult<Vec<VcfRecord>> {
    if inputs.is_empty() {
        return Err(GatkError::argument("CombineGVCFs: no inputs"));
    }
    let sample_order: Vec<String> = inputs.iter().map(|g| g.sample_name.clone()).collect();
    let mut by_contig: BTreeMap<String, Vec<(usize, &GvcfSite)>> = BTreeMap::new();
    for (si, g) in inputs.iter().enumerate() {
        for site in &g.sites {
            by_contig
                // CLONE: needed because owned HashMap entry key.
                .entry(site.contig.clone())
                .or_default()
                .push((si, site));
        }
    }

    let mut out = Vec::new();
    for (contig, sites) in by_contig {
        let clip_spans = clip.and_then(|m| m.get(&contig)).map(|v| v.as_slice());
        if clip.is_some() && clip_spans.is_none() {
            continue;
        }
        let breakpoints = collect_breakpoints(&sites, clip_spans);
        if breakpoints.len() < 2 {
            continue;
        }
        let mut contig_recs = Vec::new();
        for w in breakpoints.windows(2) {
            let loc_start = w[0];
            let loc_end = w[1].saturating_sub(1);
            if loc_end < loc_start {
                continue;
            }
            if let Some(spans) = clip_spans {
                if !spans
                    .iter()
                    .any(|s| overlaps(s.start, s.end, loc_start, loc_end))
                {
                    continue;
                }
            }

            let covering = covering_sites(&sites, loc_start);
            if covering.is_empty() {
                continue;
            }

            let mut merge_inputs = Vec::with_capacity(covering.len());
            let mut format_carry: HashMap<String, SampleData> = HashMap::new();
            for &(si, site) in &covering {
                let sample = &inputs[si].sample_name;
                // CLONE: needed because owned HashMap/BTree/HashSet key or value.
                format_carry.insert(sample.clone(), site.sample.clone());
                merge_inputs.push(site_to_merge_input(site, sample));
            }
            // Java `CombineGVCFs.endPreviousStates` walks overlapping VCs from last→first
            // before `ReferenceConfidenceVariantContextMerger.merge`, which fixes ALT order.
            merge_inputs.reverse();

            let ref_base = if merge_inputs.iter().any(|v| v.start == loc_start) {
                None
            } else {
                match reference {
                    Some(path) => Some(reference_base_at_1based(path, &contig, loc_start)?),
                    None => {
                        // Tests without FASTA: use first covering site's REF first base.
                        Some(covering[0].1.ref_allele.as_bytes()[0].to_ascii_uppercase())
                    }
                }
            };

            let merged = match merge_reference_confidence(
                &contig,
                loc_start,
                &merge_inputs,
                ref_base,
                false,
                false,
            )? {
                Some(m) => m,
                None => continue,
            };

            let rec = merge_result_to_record(&merged, loc_end, &sample_order, &format_carry)?;
            contig_recs.push(rec);
        }
        coalesce_homref_blocks(&mut contig_recs);
        out.extend(contig_recs);
    }
    Ok(out)
}

/// Load a single-sample (or first-sample) gVCF from disk.
pub fn load_gvcf(path: &Path) -> GatkResult<LoadedGvcf> {
    let mut reader = VcfReader::from_file(path)?;
    let header = reader.header().clone();
    let sample_name =
        header.samples.first().cloned().ok_or_else(|| {
            GatkError::argument(format!("gVCF has no samples: {}", path.display()))
        })?;
    let raw = reader.read_all_records()?;
    let sites = raw
        .into_iter()
        .map(|rec| record_to_site(rec, &sample_name, path))
        .collect::<GatkResult<Vec<_>>>()?;
    Ok(LoadedGvcf {
        path: path.to_path_buf(),
        sample_name,
        header,
        sites,
    })
}

/// Parse an in-memory VCF document (tests / fixtures).
#[cfg(test)]
pub fn load_gvcf_from_str(vcf_text: &str, label: &str) -> GatkResult<LoadedGvcf> {
    let dir = tempfile::tempdir().map_err(|e| GatkError::io("tempdir", e))?;
    let path = dir.path().join(format!("{label}.g.vcf"));
    std::fs::write(&path, vcf_text).map_err(|e| GatkError::io("write temp gvcf", e))?;
    let mut loaded = load_gvcf(&path)?;
    loaded.path = PathBuf::from(label);
    Ok(loaded)
}

fn resolve_clip_spans(
    dictionary: &SequenceDictionary,
    intervals: Option<&str>,
) -> GatkResult<Option<HashMap<String, Vec<Span>>>> {
    let Some(iv) = intervals else {
        return Ok(None);
    };
    let specs = gatk_core::reference::parse_intervals_cli_string(dictionary, iv)?;
    let mut map: HashMap<String, Vec<Span>> = HashMap::new();
    for spec in specs {
        let start = spec.start.unwrap_or(1);
        let end = match spec.end {
            Some(e) => e,
            None => {
                let c = dictionary.contig(&spec.contig).ok_or_else(|| {
                    GatkError::argument(format!(
                        "whole-contig -L {} needs contig length in dictionary",
                        spec.contig
                    ))
                })?;
                c.length
            }
        };
        map.entry(spec.contig)
            .or_default()
            .push(Span { start, end });
    }
    Ok(Some(map))
}

fn overlaps(a0: u64, a1: u64, b0: u64, b1: u64) -> bool {
    a0 <= b1 && b0 <= a1
}

fn record_end(rec: &VcfRecord) -> u64 {
    for info in &rec.info {
        if let InfoValue::Integer(id, vals) = info {
            if id == "END" {
                if let Some(&e) = vals.first() {
                    return e as u64;
                }
            }
        }
    }
    rec.position + rec.reference.len().saturating_sub(1) as u64
}

fn record_to_site(rec: VcfRecord, _sample_name: &str, path: &Path) -> GatkResult<GvcfSite> {
    let end = record_end(&rec);
    let sample = rec.samples.first().cloned().unwrap_or_else(|| SampleData {
        gt: None,
        gq: None,
        dp: None,
        ad: None,
        pl: None,
        other: Vec::new(),
    });
    Ok(GvcfSite {
        contig: rec.chromosome,
        start: rec.position,
        end,
        ref_allele: rec.reference,
        alts: rec.alternate,
        sample,
        source: path.display().to_string(),
    })
}

fn site_to_merge_input(site: &GvcfSite, sample: &str) -> MergeVcInput {
    let mut alleles = Vec::with_capacity(1 + site.alts.len());
    alleles.push(MergeAllele {
        bases: site.ref_allele.clone(),
        is_reference: true,
    });
    for alt in &site.alts {
        alleles.push(MergeAllele {
            // CLONE: needed because owned haplotypes for scoring call.
            bases: alt.clone(),
            is_reference: false,
        });
    }
    let pl = site.sample.pl.as_ref().map(|v| {
        v.iter()
            .map(|&x| i32::try_from(x).unwrap_or(i32::MAX))
            .collect()
    });
    let ad = site.sample.ad.as_ref().map(|v| {
        v.iter()
            .map(|&x| i32::try_from(x).unwrap_or(i32::MAX))
            .collect()
    });
    MergeVcInput {
        source: site.source.clone(),
        start: site.start,
        alleles,
        genotypes: vec![MergeGenotype {
            sample: sample.to_string(),
            pl,
            ad,
        }],
    }
}

fn collect_breakpoints(sites: &[(usize, &GvcfSite)], clip_spans: Option<&[Span]>) -> Vec<u64> {
    let mut bp = BTreeSet::new();
    for (_, s) in sites {
        bp.insert(s.start);
        bp.insert(s.end.saturating_add(1));
        // Spanning variant bodies need per-base merges (`*` remapping).
        let has_variant_alt = s.alts.iter().any(|a| a != NON_REF_ALLELE);
        if has_variant_alt && s.end > s.start {
            for p in (s.start + 1)..=s.end {
                bp.insert(p);
            }
        }
    }
    if let Some(spans) = clip_spans {
        for s in spans {
            bp.insert(s.start);
            bp.insert(s.end.saturating_add(1));
        }
    }
    bp.into_iter().collect()
}

fn covering_sites<'a>(sites: &'a [(usize, &GvcfSite)], loc: u64) -> Vec<(usize, &'a GvcfSite)> {
    let mut out = Vec::new();
    let mut seen_sample = BTreeSet::new();
    for &(si, site) in sites {
        if site.start <= loc && loc <= site.end && seen_sample.insert(si) {
            out.push((si, site));
        }
    }
    out
}

fn is_homref_nonref_only(alts: &[String]) -> bool {
    alts.len() == 1 && alts[0] == NON_REF_ALLELE
}

fn merge_result_to_record(
    merged: &crate::ref_confidence_merger::RefConfidenceMergeResult,
    interval_end: u64,
    sample_order: &[String],
    format_carry: &HashMap<String, SampleData>,
) -> GatkResult<VcfRecord> {
    let ref_allele = merged
        .alleles
        .first()
        .cloned()
        .ok_or_else(|| GatkError::argument("merged alleles empty"))?;
    let alts: Vec<String> = merged.alleles.iter().skip(1).cloned().collect();
    let n_alleles = merged.alleles.len();

    let by_name: HashMap<&str, &crate::ref_confidence_merger::MergeGenotypeOut> = merged
        .genotypes
        .iter()
        .map(|g| (g.name.as_str(), g))
        .collect();

    let mut samples = Vec::with_capacity(sample_order.len());
    for name in sample_order {
        if let Some(g) = by_name.get(name.as_str()) {
            let carry = format_carry.get(name);
            let pl_u32: Option<Vec<u32>> = g.pl.as_ref().map(|pl| {
                pl.iter()
                    .map(|&x| u32::try_from(x.max(0)).unwrap_or(u32::MAX))
                    .collect()
            });
            let ad_u32: Option<Vec<u32>> = g.ad.as_ref().map(|ad| {
                ad.iter()
                    .map(|&x| u32::try_from(x.max(0)).unwrap_or(0))
                    .collect()
            });
            let gt = pl_u32
                .as_ref()
                .and_then(|pl| best_diploid_gt_from_pl_u32(pl.as_slice(), n_alleles));
            // CLONE: needed because multi-owner or ownership transfer into new structure.
            let other = carry.map(|c| c.other.clone()).unwrap_or_default();
            // Preserve MIN_DP and any non-standard FORMAT keys from the covering record.
            samples.push(SampleData {
                gt,
                gq: carry.and_then(|c| c.gq),
                dp: carry.and_then(|c| c.dp),
                ad: ad_u32,
                pl: pl_u32,
                other,
            });
        } else {
            samples.push(SampleData {
                gt: None,
                gq: None,
                dp: None,
                ad: None,
                pl: None,
                other: Vec::new(),
            });
        }
    }

    let mut info = Vec::new();
    let block = is_homref_nonref_only(&alts) && ref_allele.len() == 1;
    if block && interval_end > merged.pos {
        info.push(InfoValue::Integer(
            "END".to_string(),
            vec![interval_end as i32],
        ));
    }

    Ok(VcfRecord {
        chromosome: merged.contig.clone(),
        position: merged.pos,
        id: ".".to_string(),
        reference: ref_allele,
        alternate: alts,
        quality: None,
        filter: vec![".".to_string()],
        info,
        format: vec![
            "GT".to_string(),
            "GQ".to_string(),
            "DP".to_string(),
            "AD".to_string(),
            "PL".to_string(),
            "MIN_DP".to_string(),
        ],
        samples,
    })
}

fn best_diploid_gt_from_pl_u32(pl: &[u32], n_alleles: usize) -> Option<Genotype> {
    if n_alleles == 0 || pl.is_empty() {
        return None;
    }
    let mut pairs = Vec::new();
    for j in 0..n_alleles {
        for i in 0..=j {
            pairs.push((i, j));
        }
    }
    if pl.len() != pairs.len() {
        return None;
    }
    let (best_idx, _) = pl.iter().enumerate().min_by_key(|(_, &v)| v)?;
    let (a, b) = pairs[best_idx];
    Some(Genotype {
        alleles: vec![a as i32, b as i32],
        phased: false,
    })
}

/// Merge adjacent hom-ref `<NON_REF>` blocks with identical per-sample FORMAT payloads.
fn coalesce_homref_blocks(recs: &mut Vec<VcfRecord>) {
    if recs.len() < 2 {
        return;
    }
    let mut out: Vec<VcfRecord> = Vec::with_capacity(recs.len());
    for rec in recs.drain(..) {
        let can_merge = out.last().is_some_and(|prev| {
            is_homref_nonref_only(&prev.alternate)
                && is_homref_nonref_only(&rec.alternate)
                && prev.chromosome == rec.chromosome
                && prev.reference == rec.reference
                && samples_pl_ad_eq(&prev.samples, &rec.samples)
                && record_end(prev) + 1 == rec.position
        });
        if can_merge {
            let new_end = record_end(&rec);
            if let Some(prev) = out.last_mut() {
                set_end(prev, new_end);
            }
        } else {
            out.push(rec);
        }
    }
    *recs = out;
}

fn samples_pl_ad_eq(a: &[SampleData], b: &[SampleData]) -> bool {
    if a.len() != b.len() {
        return false;
    }
    a.iter().zip(b.iter()).all(|(x, y)| {
        x.pl == y.pl && x.ad == y.ad && x.dp == y.dp && x.gq == y.gq && x.other == y.other
    })
}

fn set_end(rec: &mut VcfRecord, end: u64) {
    rec.info
        .retain(|v| !matches!(v, InfoValue::Integer(id, _) if id == "END"));
    if end > rec.position {
        rec.info
            .push(InfoValue::Integer("END".to_string(), vec![end as i32]));
    }
}

fn build_combined_header(inputs: &[LoadedGvcf], reference: Option<&Path>) -> GatkResult<VcfHeader> {
    let mut header = VcfHeader::default();
    header.file_format = "VCFv4.2".to_string();
    header.source = Some("gatk-rs CombineGVCFs".to_string());
    if let Some(r) = reference {
        header.reference = Some(r.display().to_string());
    }

    let mut contig_ids = BTreeSet::new();
    for g in inputs {
        for c in &g.header.contigs {
            if c.id.is_empty() {
                continue;
            }
            // CLONE: needed because owned HashMap/BTree/HashSet key or value.
            if contig_ids.insert(c.id.clone()) {
                // CLONE: needed because owned element into collection.
                header.contigs.push(c.clone());
            }
        }
    }
    if header.contigs.is_empty() {
        let mut seen = BTreeSet::new();
        for g in inputs {
            for s in &g.sites {
                // CLONE: needed because owned HashMap/BTree/HashSet key or value.
                if seen.insert(s.contig.clone()) {
                    header.contigs.push(Contig {
                        // CLONE: needed because VCF header contig id is owned separately from `seen`.
                        id: s.contig.clone(),
                        length: None,
                        md5: None,
                        assembly: None,
                        species: None,
                        uri: None,
                    });
                }
            }
        }
    }

    header.info_fields = vec![
        InfoField {
            id: "END".to_string(),
            number: "1".to_string(),
            type_field: "Integer".to_string(),
            description: "Stop position of the interval".to_string(),
            source: None,
            version: None,
        },
        InfoField {
            id: "DP".to_string(),
            number: "1".to_string(),
            type_field: "Integer".to_string(),
            description: "Approximate read depth; some reads may have been filtered".to_string(),
            source: None,
            version: None,
        },
    ];
    header.format_fields = vec![
        FormatField {
            id: "GT".to_string(),
            number: "1".to_string(),
            type_field: "String".to_string(),
            description: "Genotype".to_string(),
        },
        FormatField {
            id: "GQ".to_string(),
            number: "1".to_string(),
            type_field: "Integer".to_string(),
            description: "Genotype Quality".to_string(),
        },
        FormatField {
            id: "DP".to_string(),
            number: "1".to_string(),
            type_field: "Integer".to_string(),
            description:
                "Approximate read depth (reads with MQ=255 or with bad mates are filtered)"
                    .to_string(),
        },
        FormatField {
            id: "AD".to_string(),
            number: "R".to_string(),
            type_field: "Integer".to_string(),
            description: "Allelic depths for the ref and alt alleles in the order listed"
                .to_string(),
        },
        FormatField {
            id: "PL".to_string(),
            number: "G".to_string(),
            type_field: "Integer".to_string(),
            description: "Normalized, Phred-scaled likelihoods for genotypes".to_string(),
        },
        FormatField {
            id: "MIN_DP".to_string(),
            number: "1".to_string(),
            type_field: "Integer".to_string(),
            description: "Minimum DP observed within the GVCF block".to_string(),
        },
    ];
    header.samples = inputs.iter().map(|g| g.sample_name.clone()).collect();
    Ok(header)
}

// --------------------------------------------------------------------------
// Tests
// --------------------------------------------------------------------------

#[cfg(test)]
#[allow(clippy::unwrap_used, clippy::expect_used)]
mod tests {
    use super::*;
    use crate::ref_confidence_merger::SPAN_DEL_ALLELE;

    fn mini_header(sample: &str, contig: &str, contig_len: u64) -> String {
        format!(
            "##fileformat=VCFv4.2\n\
             ##contig=<ID={contig},length={contig_len}>\n\
             ##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop\">\n\
             ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n\
             ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"GQ\">\n\
             ##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"DP\">\n\
             ##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"AD\">\n\
             ##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"PL\">\n\
             ##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=\"MIN_DP\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
        )
    }

    fn load_two(a: &str, b: &str) -> (LoadedGvcf, LoadedGvcf) {
        (
            load_gvcf_from_str(a, "sA").unwrap(),
            load_gvcf_from_str(b, "sB").unwrap(),
        )
    }

    fn pl_of<'a>(rec: &'a VcfRecord, sample_idx: usize) -> &'a [u32] {
        rec.samples[sample_idx].pl.as_deref().unwrap_or(&[])
    }

    #[test]
    fn t01_two_identical_ref_blocks_merge() {
        let h = mini_header("S1", "chr1", 200);
        let a = format!(
            "{h}chr1\t10\t.\tA\t<NON_REF>\t.\t.\tEND=20\tGT:GQ:DP:MIN_DP:PL\t0/0:99:30:25:0,90,900\n"
        );
        let h2 = mini_header("S2", "chr1", 200);
        let b = format!(
            "{h2}chr1\t10\t.\tA\t<NON_REF>\t.\t.\tEND=20\tGT:GQ:DP:MIN_DP:PL\t0/0:80:20:18:0,60,600\n"
        );
        let (la, lb) = load_two(&a, &b);
        let out = combine_loaded_gvcfs(&[la, lb], None, None).unwrap();
        assert!(!out.is_empty());
        assert_eq!(out[0].position, 10);
        assert_eq!(record_end(&out[0]), 20);
        assert_eq!(out[0].alternate, vec![NON_REF_ALLELE]);
        assert_eq!(out[0].samples.len(), 2);
        assert_eq!(pl_of(&out[0], 0), &[0, 90, 900]);
        assert_eq!(pl_of(&out[0], 1), &[0, 60, 600]);
    }

    #[test]
    fn t02_overlapping_ref_blocks_different_lengths() {
        let a = format!(
            "{}chr1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=10\tGT:DP:PL\t0/0:10:0,30,300\n",
            mini_header("A", "chr1", 100)
        );
        let b = format!(
            "{}chr1\t5\t.\tC\t<NON_REF>\t.\t.\tEND=15\tGT:DP:PL\t0/0:12:0,36,360\n",
            mini_header("B", "chr1", 100)
        );
        let (la, lb) = load_two(&a, &b);
        let out = combine_loaded_gvcfs(&[la, lb], None, None).unwrap();
        // Breakpoints at 1,5,11,16 → intervals covering both overlap and exclusive tails.
        assert!(out.iter().any(|r| r.position == 1));
        assert!(out.iter().any(|r| r.position == 5));
        let at5 = out.iter().find(|r| r.position == 5).unwrap();
        assert_eq!(at5.samples.len(), 2);
        // Sample A still covered (spanning), sample B starts here.
        assert!(at5.samples[0].pl.is_some());
        assert!(at5.samples[1].pl.is_some());
    }

    #[test]
    fn t03_sample_without_coverage_is_nocall() {
        let a = format!(
            "{}chr1\t10\t.\tA\t<NON_REF>\t.\t.\tEND=12\tGT:DP:PL\t0/0:5:0,15,150\n",
            mini_header("A", "chr1", 100)
        );
        let b = format!(
            "{}chr1\t50\t.\tT\t<NON_REF>\t.\t.\tEND=52\tGT:DP:PL\t0/0:5:0,15,150\n",
            mini_header("B", "chr1", 100)
        );
        let (la, lb) = load_two(&a, &b);
        let out = combine_loaded_gvcfs(&[la, lb], None, None).unwrap();
        let at10 = out.iter().find(|r| r.position == 10).unwrap();
        assert!(at10.samples[0].pl.is_some());
        assert!(at10.samples[1].pl.is_none());
        assert!(at10.samples[1].gt.is_none());
    }

    #[test]
    fn t04_different_snp_alts_unify_alleles() {
        // Exact mini-cohort PL failure case (parity log 20260723T221851Z, chr1:10):
        // SAMPLE1 A→G, SAMPLE2 A→T. Java `endPreviousStates` last→first ⇒ ALT order T,G,<NON_REF>.
        // Hand-calculated remap (see ref_confidence_merger::pl_remap_tests::pl01):
        // SAMPLE1 PL → 100,100,100,0,100,100,100,100,100,100 (zero at 0/2 = REF/G)
        // SAMPLE2 PL → 90,0,90,90,90,90,90,90,90,90 (zero at 0/1 = REF/T)
        let a = format!(
            "{}chr1\t10\t.\tA\tG,<NON_REF>\t.\t.\t.\tGT:AD:DP:GQ:PL\t0/1:10,10,0:20:99:100,0,100,100,100,100\n",
            mini_header("SAMPLE1", "chr1", 200)
        );
        let b = format!(
            "{}chr1\t10\t.\tA\tT,<NON_REF>\t.\t.\t.\tGT:AD:DP:GQ:PL\t0/1:8,8,0:16:90:90,0,90,90,90,90\n",
            mini_header("SAMPLE2", "chr1", 200)
        );
        let (la, lb) = load_two(&a, &b);
        let out = combine_loaded_gvcfs(&[la, lb], None, None).unwrap();
        let site = out.iter().find(|r| r.position == 10).unwrap();
        assert_eq!(
            site.alternate,
            vec!["T".to_string(), "G".to_string(), NON_REF_ALLELE.to_string()]
        );
        assert_eq!(
            pl_of(site, 0),
            &[100, 100, 100, 0, 100, 100, 100, 100, 100, 100]
        );
        assert_eq!(pl_of(site, 1), &[90, 0, 90, 90, 90, 90, 90, 90, 90, 90]);
    }

    #[test]
    fn t05_indel_spanning_into_ref_block() {
        // A: deletion ACG→A at 10 (spans 10-12). B: ref block covering 10-20.
        let a = format!(
            "{}chr1\t10\t.\tACG\tA,<NON_REF>\t.\t.\t.\tGT:AD:DP:PL\t0/1:5,5,0:10:50,0,50,50,50,50\n",
            mini_header("A", "chr1", 100)
        );
        let b = format!(
            "{}chr1\t10\t.\tA\t<NON_REF>\t.\t.\tEND=20\tGT:DP:PL\t0/0:8:0,24,240\n",
            mini_header("B", "chr1", 100)
        );
        let (la, lb) = load_two(&a, &b);
        let out = combine_loaded_gvcfs(&[la, lb], None, None).unwrap();
        let at10 = out.iter().find(|r| r.position == 10).unwrap();
        assert!(at10
            .alternate
            .iter()
            .any(|a| a == "A" || a.starts_with('A')));
        // Interior of deletion should appear as subsequent breakpoints.
        assert!(out.iter().any(|r| r.position == 11 || r.position == 12));
    }

    #[test]
    fn t06_spanning_deletion_emits_star_allele() {
        let a = format!(
            "{}chr1\t10\t.\tACG\tA,<NON_REF>\t.\t.\t.\tGT:AD:PL\t0/1:4,4,0:40,0,40,40,40,40\n",
            mini_header("A", "chr1", 100)
        );
        let b = format!(
            "{}chr1\t11\t.\tC\tT,<NON_REF>\t.\t.\t.\tGT:AD:PL\t0/1:3,3,0:30,0,30,30,30,30\n",
            mini_header("B", "chr1", 100)
        );
        let (la, lb) = load_two(&a, &b);
        let out = combine_loaded_gvcfs(&[la, lb], None, None).unwrap();
        let at11 = out.iter().find(|r| r.position == 11).unwrap();
        assert!(
            at11.alternate.iter().any(|a| a == SPAN_DEL_ALLELE),
            "expected * from spanning indel, got {:?}",
            at11.alternate
        );
    }

    #[test]
    fn t07_preserves_min_dp_in_other() {
        let a = format!(
            "{}chr1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=5\tGT:DP:MIN_DP:PL\t0/0:40:22:0,99,999\n",
            mini_header("A", "chr1", 50)
        );
        let b = format!(
            "{}chr1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=5\tGT:DP:MIN_DP:PL\t0/0:10:7:0,30,300\n",
            mini_header("B", "chr1", 50)
        );
        let (la, lb) = load_two(&a, &b);
        let out = combine_loaded_gvcfs(&[la, lb], None, None).unwrap();
        let min_dp = |s: &SampleData| -> Option<String> {
            s.other
                .iter()
                .find(|(k, _)| k == "MIN_DP")
                .map(|(_, v)| v.clone())
        };
        assert_eq!(min_dp(&out[0].samples[0]).as_deref(), Some("22"));
        assert_eq!(min_dp(&out[0].samples[1]).as_deref(), Some("7"));
    }

    #[test]
    fn t08_three_samples_union() {
        let mk = |sample: &str, dp: u32| {
            format!(
                "{}chr1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=3\tGT:DP:PL\t0/0:{dp}:0,10,100\n",
                mini_header(sample, "chr1", 20)
            )
        };
        let inputs = [
            load_gvcf_from_str(&mk("X", 1), "X").unwrap(),
            load_gvcf_from_str(&mk("Y", 2), "Y").unwrap(),
            load_gvcf_from_str(&mk("Z", 3), "Z").unwrap(),
        ];
        let out = combine_loaded_gvcfs(&inputs, None, None).unwrap();
        assert_eq!(out[0].samples.len(), 3);
        assert_eq!(out[0].samples[0].dp, Some(1));
        assert_eq!(out[0].samples[2].dp, Some(3));
    }

    #[test]
    fn t09_clip_interval_excludes_outside() {
        let a = format!(
            "{}chr1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=100\tGT:DP:PL\t0/0:5:0,15,150\n",
            mini_header("A", "chr1", 200)
        );
        let b = format!(
            "{}chr1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=100\tGT:DP:PL\t0/0:5:0,15,150\n",
            mini_header("B", "chr1", 200)
        );
        let (la, lb) = load_two(&a, &b);
        let mut clip = HashMap::new();
        clip.insert("chr1".to_string(), vec![Span { start: 50, end: 60 }]);
        let out = combine_loaded_gvcfs(&[la, lb], None, Some(&clip)).unwrap();
        assert!(out.iter().all(|r| r.position >= 50 && record_end(r) <= 60));
    }

    #[test]
    fn t10_coalesce_adjacent_identical_homref() {
        // Two abutting blocks with same FORMAT → one coalesced END.
        let a = format!(
            "{}chr1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=5\tGT:DP:PL\t0/0:5:0,15,150\n\
             chr1\t6\t.\tC\t<NON_REF>\t.\t.\tEND=10\tGT:DP:PL\t0/0:5:0,15,150\n",
            mini_header("A", "chr1", 50)
        );
        let b = format!(
            "{}chr1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=5\tGT:DP:PL\t0/0:5:0,15,150\n\
             chr1\t6\t.\tC\t<NON_REF>\t.\t.\tEND=10\tGT:DP:PL\t0/0:5:0,15,150\n",
            mini_header("B", "chr1", 50)
        );
        let (la, lb) = load_two(&a, &b);
        let out = combine_loaded_gvcfs(&[la, lb], None, None).unwrap();
        // After coalesce, prefer fewer hom-ref rows than raw breakpoints.
        let homref = out
            .iter()
            .filter(|r| is_homref_nonref_only(&r.alternate))
            .count();
        assert!(homref >= 1);
        assert!(out.iter().any(|r| record_end(r) >= 5));
    }

    #[test]
    fn t11_allele_remap_maps_missing_alt_to_nonref_pl() {
        let a = format!(
            "{}chr1\t1\t.\tA\tG,<NON_REF>\t.\t.\t.\tGT:AD:PL\t0/0:10,0,0:0,100,1000,100,1000,1000\n",
            mini_header("A", "chr1", 10)
        );
        let b = format!(
            "{}chr1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=1\tGT:DP:PL\t0/0:10:0,30,300\n",
            mini_header("B", "chr1", 10)
        );
        let (la, lb) = load_two(&a, &b);
        let out = combine_loaded_gvcfs(&[la, lb], None, None).unwrap();
        let site = &out[0];
        assert!(site.alternate.contains(&"G".to_string()));
        // B had 3 PL values; remapped onto REF,G,<NON_REF> → 6 PLs.
        assert_eq!(pl_of(site, 1).len(), 6);
    }

    #[test]
    fn t12_empty_inputs_error() {
        let err = combine_loaded_gvcfs(&[], None, None).unwrap_err();
        assert!(format!("{err}").contains("no inputs"));
    }

    #[test]
    fn t13_record_end_from_info_and_ref_len() {
        let with_end = format!(
            "{}chr1\t10\t.\tA\t<NON_REF>\t.\t.\tEND=42\tGT:PL\t0/0:0,1,2\n",
            mini_header("A", "chr1", 100)
        );
        let indel = format!(
            "{}chr1\t10\t.\tACGT\tA,<NON_REF>\t.\t.\t.\tGT:PL\t0/1:10,0,10,10,10,10\n",
            mini_header("B", "chr1", 100)
        );
        let la = load_gvcf_from_str(&with_end, "A").unwrap();
        let lb = load_gvcf_from_str(&indel, "B").unwrap();
        assert_eq!(la.sites[0].end, 42);
        assert_eq!(lb.sites[0].end, 13); // 10 + 4 - 1
    }

    #[test]
    fn t14_gq_dp_carried_from_input() {
        let a = format!(
            "{}chr1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=2\tGT:GQ:DP:PL\t0/0:55:17:0,20,200\n",
            mini_header("A", "chr1", 10)
        );
        let b = format!(
            "{}chr1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=2\tGT:GQ:DP:PL\t0/0:40:9:0,12,120\n",
            mini_header("B", "chr1", 10)
        );
        let (la, lb) = load_two(&a, &b);
        let out = combine_loaded_gvcfs(&[la, lb], None, None).unwrap();
        assert_eq!(out[0].samples[0].gq, Some(55.0));
        assert_eq!(out[0].samples[0].dp, Some(17));
        assert_eq!(out[0].samples[1].gq, Some(40.0));
        assert_eq!(out[0].samples[1].dp, Some(9));
    }

    #[test]
    fn t15_variant_site_has_no_end_when_alts_present() {
        let a = format!(
            "{}chr1\t10\t.\tA\tG,<NON_REF>\t.\t.\t.\tGT:AD:PL\t0/1:5,5,0:50,0,50,50,50,50\n",
            mini_header("A", "chr1", 50)
        );
        let b = format!(
            "{}chr1\t10\t.\tA\tG,<NON_REF>\t.\t.\t.\tGT:AD:PL\t0/1:4,4,0:40,0,40,40,40,40\n",
            mini_header("B", "chr1", 50)
        );
        let (la, lb) = load_two(&a, &b);
        let out = combine_loaded_gvcfs(&[la, lb], None, None).unwrap();
        let site = out.iter().find(|r| r.position == 10).unwrap();
        assert!(site.alternate.contains(&"G".to_string()));
        assert!(!site
            .info
            .iter()
            .any(|i| matches!(i, InfoValue::Integer(id, _) if id == "END")));
    }

    #[test]
    fn t16_ad_remapped_length_matches_alleles() {
        let a = format!(
            "{}chr1\t1\t.\tA\tC,<NON_REF>\t.\t.\t.\tGT:AD:PL\t0/1:7,3,0:30,0,90,90,90,90\n",
            mini_header("A", "chr1", 10)
        );
        let b = format!(
            "{}chr1\t1\t.\tA\tG,<NON_REF>\t.\t.\t.\tGT:AD:PL\t0/1:6,4,0:40,0,80,80,80,80\n",
            mini_header("B", "chr1", 10)
        );
        let (la, lb) = load_two(&a, &b);
        let out = combine_loaded_gvcfs(&[la, lb], None, None).unwrap();
        let n_alleles = 1 + out[0].alternate.len();
        assert_eq!(out[0].samples[0].ad.as_ref().unwrap().len(), n_alleles);
        assert_eq!(out[0].samples[1].ad.as_ref().unwrap().len(), n_alleles);
    }

    #[test]
    fn t17_multi_contig_keeps_order() {
        let a = format!(
            "{}chr1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=2\tGT:PL\t0/0:0,1,2\n\
             chr2\t1\t.\tT\t<NON_REF>\t.\t.\tEND=2\tGT:PL\t0/0:0,1,2\n",
            mini_header("A", "chr1", 10).replace(
                "##contig=<ID=chr1,length=10>",
                "##contig=<ID=chr1,length=10>\n##contig=<ID=chr2,length=10>"
            )
        );
        let b = format!(
            "{}chr1\t1\t.\tA\t<NON_REF>\t.\t.\tEND=2\tGT:PL\t0/0:0,1,2\n\
             chr2\t1\t.\tT\t<NON_REF>\t.\t.\tEND=2\tGT:PL\t0/0:0,1,2\n",
            mini_header("B", "chr1", 10).replace(
                "##contig=<ID=chr1,length=10>",
                "##contig=<ID=chr1,length=10>\n##contig=<ID=chr2,length=10>"
            )
        );
        let (la, lb) = load_two(&a, &b);
        let out = combine_loaded_gvcfs(&[la, lb], None, None).unwrap();
        let contigs: Vec<_> = out.iter().map(|r| r.chromosome.as_str()).collect();
        assert!(contigs.contains(&"chr1"));
        assert!(contigs.contains(&"chr2"));
        let first_chr2 = contigs.iter().position(|c| *c == "chr2").unwrap();
        let last_chr1 = contigs.iter().rposition(|c| *c == "chr1").unwrap();
        assert!(last_chr1 < first_chr2);
    }
}
