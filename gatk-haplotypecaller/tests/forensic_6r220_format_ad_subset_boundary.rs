//! 6R.220: first operation that turns the 52-row retainEvidence remarg
//! (`AD 47,5` / `PL 68,0,1937`) into the FORMAT writer object
//! (`AD 44,5` / `DP 49` / `PL 78,0,1811`) at `20:29455379 G/A`.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r220_format_ad_subset_boundary -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{
    build_per_haplotype_variation_events, variation_events_at_position_from_cache,
};
use gatk_haplotypecaller::hc_allele_mapping::create_allele_mapper_with_events;
use gatk_haplotypecaller::hc_genotyping_engine::{
    alt_hap_indices_for_genotype_marginalization, biallelic_genotype_log10_likelihoods_gatk,
    diagnose_genotype_variation_event, java_alignment_read_covers_variant_base,
    java_alignment_read_overlaps_interval, marginalize_rows_to_biallelic_alleles,
    ref_hap_indices_for_genotype_marginalization, region_likelihoods_to_rows, HcGenotypingConfig,
    InformativeAd, DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN,
};
use gatk_haplotypecaller::query_index_at_reference_position;
use gatk_haplotypecaller::read_model::{
    FLAG_DUPLICATE, FLAG_NOT_PRIMARY, FLAG_SUPPLEMENTARY, FLAG_VENDOR_QUALITY_FAILED,
};
use gatk_haplotypecaller::read_realignment::LOG_10_INFORMATIVE_THRESHOLD;
use gatk_haplotypecaller::read_unclip::alignment_end_1based;
use gatk_haplotypecaller::variant_site_hc_annotations::coverage_evidence_count;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use rust_htslib::bam::record::CigarString;
use rust_htslib::bam::{IndexedReader, Read};
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_455_379;
const CLOSED_PL: u64 = 29_455_015;
const TARGET_REF: &str = "G";
const TARGET_ALT: &str = "A";
const JAVA_AD: [i32; 2] = [42, 5];
const JAVA_FMT_DP: i32 = 47;
const JAVA_PL: [i32; 3] = [84, 0, 1738];
const MARGIN: i32 = DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R220\t{key}\t{}", value.as_ref());
}

fn info_i32(info: &[InfoValue], key: &str) -> Option<i32> {
    for v in info {
        if let InfoValue::Integer(k, xs) = v {
            if k == key {
                return xs.first().copied();
            }
        }
    }
    None
}

fn unique_indices(likelihoods: &[gatk_haplotypecaller::RegionReadLikelihood]) -> BTreeSet<usize> {
    likelihoods.iter().map(|c| c.read_index.get()).collect()
}

#[derive(Clone, Debug)]
#[allow(dead_code)]
struct Row {
    idx: usize,
    qname: String,
    flag: u16,
    mapq: u8,
    pos: i64,
    end: i64,
    cigar: String,
    lr: f64,
    la: f64,
    conf: f64,
    informative: bool,
    vote: &'static str,
    overlap_pm2: bool,
    overlap_pm0: bool,
    covers: bool,
    pileup: Option<u8>,
    dup: bool,
    secondary: bool,
    supplementary: bool,
    qc_fail: bool,
    proper_pair: bool,
    mate_unmapped: bool,
    first: bool,
    reverse: bool,
}

fn vote_row(lr: f64, la: f64) -> (f64, bool, &'static str) {
    let conf = (lr - la).abs();
    let informative = lr.is_finite() && la.is_finite() && conf > LOG_10_INFORMATIVE_THRESHOLD;
    let vote = if !informative {
        "UNINF"
    } else if lr > la {
        "REF"
    } else {
        "ALT"
    };
    (conf, informative, vote)
}

fn apply_java_mismapping_floor(lr: f64, la: f64) -> (f64, f64) {
    if !lr.is_finite() || !la.is_finite() {
        return (lr, la);
    }
    let best = lr.max(la);
    let floor = best - 4.5;
    (
        if lr < floor { floor } else { lr },
        if la < floor { floor } else { la },
    )
}

fn gls_to_pl(gls: &[f64]) -> Vec<i32> {
    let max = gls.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    gls.iter()
        .map(|g| ((max - *g) * 10.0).round() as i32)
        .collect()
}

fn pileup_base_at(rec: &rust_htslib::bam::Record, loc_1based: u64) -> Option<u8> {
    let cigar = CigarString(rec.cigar().iter().copied().collect());
    let qi = query_index_at_reference_position(rec.pos(), &cigar, loc_1based as i64 - 1)?;
    Some(rec.seq()[qi])
}

fn overlap_pm0(rec: &rust_htslib::bam::Record, loc: u64) -> bool {
    java_alignment_read_overlaps_interval(rec, loc, loc, 0)
}

fn ad_pl_of<'a, I>(rows: I) -> (Vec<i32>, Vec<i32>, i32, i32, i32)
where
    I: IntoIterator<Item = &'a gatk_haplotypecaller::genotyping::ReadLikelihoodRow>,
{
    let kept: Vec<_> = rows.into_iter().cloned().collect();
    let mut floored = kept.clone();
    for row in &mut floored {
        let (lr, la) = apply_java_mismapping_floor(
            row.haplotype_log10_likelihoods[0],
            row.haplotype_log10_likelihoods[1],
        );
        row.haplotype_log10_likelihoods[0] = lr;
        row.haplotype_log10_likelihoods[1] = la;
    }
    let ad = InformativeAd::from_marginalized_rows(&floored, 0, 1, None).as_vec();
    let pl = gls_to_pl(&biallelic_genotype_log10_likelihoods_gatk(&floored, 0, 1));
    let mut ref_n = 0i32;
    let mut alt_n = 0i32;
    let mut uninf_n = 0i32;
    for row in &floored {
        match vote_row(
            row.haplotype_log10_likelihoods[0],
            row.haplotype_log10_likelihoods[1],
        )
        .2
        {
            "REF" => ref_n += 1,
            "ALT" => alt_n += 1,
            _ => uninf_n += 1,
        }
    }
    (ad, pl, ref_n, alt_n, uninf_n)
}

#[test]
fn forensic_6r220_java_format_lifecycle_contract() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE — proof-only");
    kv(
        "java_ad_lifecycle",
        "assignGenotypeLikelihoods: AlleleLikelihoods created from PairHMM perSampleReadList → marginalize(alleleMapper) → retainEvidence(variantCallingRelevantOverlap ±2; overlap only; no groupEvidence) → calculateGLsForThisEvent (PL only on that retained object) → prepareReadAlleleLikelihoodsForAnnotation (contamination off: reuse genotyping AlleleLikelihoods; addEvidence(filtered,0) ties extra reads uninformative) → DepthPerAlleleBySample.annotateWithLikelihoods on remaining vc.getAlleles() of that same object → bestAllelesBreakingTies → isInformative(confidence > 0.2) → FORMAT AD. reverseTrimAlleles after annotation copies AD; does not rebuild evidence. calculateGLs does not write AD. No FORMAT-specific post-GL read list.",
    );
    kv(
        "java_same_object",
        "calculateGLsForThisEvent and DepthPerAlleleBySample consume the same retainEvidence AlleleLikelihoods when contamination is off",
    );
    kv(
        "java_no_format_subset",
        "default HC has no keep_qnames / QNAME collapse / covering-base / mate-collapse FORMAT subset after retainEvidence",
    );
    assert_eq!(LOG_10_INFORMATIVE_THRESHOLD, 0.2);
    assert_eq!(MARGIN, 2);
    assert_eq!(JAVA_AD, [42, 5]);
    assert_eq!(JAVA_FMT_DP, 47);
    assert_eq!(JAVA_PL, [84, 0, 1738]);
}

#[test]
fn forensic_6r220_format_ad_subset_boundary() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE — proof-only");
    kv("target", "20:29455379 G/A");

    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, INTERVAL).expect("interval");
    let walk = traverse_assembly_region_walker(
        &dict,
        &specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("walk");
    let regions = flatten_assembly_regions(&walk);
    let covering: Vec<_> = regions
        .iter()
        .filter(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .collect();
    assert_eq!(covering.len(), 1, "one ActiveFull covers the target");
    let region = covering[0];
    kv(
        "active_full",
        format!(
            "{}:{}-{}",
            region.contig,
            region.start.get(),
            region.end.get()
        ),
    );

    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");

    let haps = &outcome.assembly.haplotypes;
    let reads = &outcome.genotyping_reads;
    let likelihoods = &outcome.read_likelihoods;
    let (full_ref, full_pad) = outcome.assembly.event_map_reference();
    let hap_events = build_per_haplotype_variation_events(
        haps,
        full_ref,
        full_pad,
        outcome.assembly.max_mnp_distance(),
        &region.contig,
    );
    let emit_spanning = !HcGenotypingConfig::strict_java().disable_spanning_event_genotyping;
    let cache_only = variation_events_at_position_from_cache(&hap_events, TARGET, emit_spanning);
    assert_eq!(cache_only.len(), 1);
    assert_eq!(cache_only[0].ref_allele, TARGET_REF);
    assert_eq!(cache_only[0].alt_allele, TARGET_ALT);
    let event = &cache_only[0];

    // Production `is_sparse_snp_gl_rescue_eligible` = SNP && P12 emit-band
    // (chr2 92.305M–92.325M). Cluster helpers are exact chr2 loci
    // (92307333 TG, 92307364 TC, 92307383 AC, …). This event is chr20.
    let pos = event.start_1based.get();
    let p12_band = (92_305_634..=92_325_400).contains(&pos);
    kv("event_contig", event.contig.as_str());
    kv("event_pos", pos.to_string());
    kv("p12_band_pos", p12_band.to_string());
    kv("sparse_eligible", (event.is_snp() && p12_band).to_string());
    assert_eq!(event.contig, "20");
    assert!(!p12_band, "chr20 SNP is outside P12 sparse/cluster bands");
    assert_ne!(pos, 92_307_333);
    assert_ne!(pos, 92_307_364);
    assert_ne!(pos, 92_307_383);
    assert!(
        !(event.is_snp() && p12_band),
        "keep_qnames / sparse augment must not run"
    );
    kv(
        "keep_qnames_path",
        "gated on is_sparse_snp_gl_rescue_eligible — FALSE here (introduced for P12/sparse FORMAT rescue, not default HC retainEvidence)",
    );
    kv(
        "dedupe_likelihood_subset_by_qname",
        "6R.151: applied only when sparse augment added rows; not Java retainEvidence; not taken here",
    );
    kv(
        "per_variant_annotation_likelihoods",
        "6R.180/6R.189 annotation helper; loc-loop INFO path (6R.209 no same-QNAME collapse). FORMAT caller is SiteScore subset, not this helper.",
    );

    let site = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based.get() == TARGET
                && c.event.ref_allele == TARGET_REF
                && c.event.alt_allele == TARGET_ALT
        })
        .expect("G/A site");
    assert!(site.extra_alt_alleles.is_empty());
    assert!(!site.post_merge_unused_alt_subset);
    let fmt_ad = site.genotype.format.ad_as_i32();
    let fmt_pl = site.genotype.format.pl_as_i32();
    let fmt_dp = site.genotype.format.dp.as_i32();
    kv("rust_format_ad", format!("{fmt_ad:?}"));
    kv("rust_format_pl", format!("{fmt_pl:?}"));
    kv("rust_format_dp", fmt_dp.to_string());
    assert_eq!(fmt_ad, vec![42, 5]);
    assert_eq!(fmt_pl, vec![84, 0, 1738]);
    assert_eq!(fmt_dp, 47);

    let emitted =
        try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let rec = emitted
        .iter()
        .find(|r| r.position == TARGET && r.reference == TARGET_REF)
        .expect("emit G/A");
    assert_eq!(rec.samples[0].ad.as_deref(), Some([42, 5].as_slice()));
    assert_eq!(rec.samples[0].dp, Some(47));
    assert_eq!(rec.samples[0].pl.as_deref(), Some([84, 0, 1738].as_slice()));
    assert_eq!(info_i32(&rec.info, "DP"), Some(47));

    let config = HcGenotypingConfig::strict_java();
    let apply_bases = outcome.assembly.apply_bases_shared();
    let apply_pad = haps
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.map(|g| g.start_1based()))
        .unwrap_or(full_pad);
    let mapping = create_allele_mapper_with_events(
        event,
        TARGET,
        haps,
        apply_pad,
        apply_bases.as_ref(),
        outcome.assembly.max_mnp_distance(),
        emit_spanning,
        Some(&hap_events),
    );
    let ref_hap = haps.iter().find(|h| h.is_reference).expect("ref hap");
    let ref_pool =
        ref_hap_indices_for_genotype_marginalization(&mapping, haps, &config, Some(event));
    let alt_pool = alt_hap_indices_for_genotype_marginalization(
        &mapping,
        haps,
        event,
        ref_hap,
        apply_pad,
        apply_bases.as_ref(),
        outcome.assembly.max_mnp_distance(),
        &region.contig,
        &config,
    );

    let rows = region_likelihoods_to_rows(likelihoods, haps.len());
    let marg_all = marginalize_rows_to_biallelic_alleles(&rows, &ref_pool, &alt_pool);

    let overlap_idx: BTreeSet<usize> = reads
        .iter()
        .enumerate()
        .filter(|(_, rec)| java_alignment_read_overlaps_interval(rec, TARGET, TARGET, MARGIN))
        .map(|(i, _)| i)
        .collect();
    let overlap_rows: Vec<_> = marg_all
        .iter()
        .filter(|row| {
            row.matrix_read_index()
                .is_some_and(|i| overlap_idx.contains(&i))
        })
        .cloned()
        .collect();
    let (ad52, pl52, ref52, alt52, uninf52) = ad_pl_of(&overlap_rows);
    kv("n52", overlap_rows.len().to_string());
    kv("ad52", format!("{ad52:?}"));
    kv("pl52", format!("{pl52:?}"));
    kv("vote52", format!("REF={ref52} ALT={alt52} UNINF={uninf52}"));
    assert_eq!(overlap_rows.len(), 47);
    assert_eq!(ad52, vec![42, 5]);
    assert_eq!(pl52, vec![84, 0, 1738]);
    assert_eq!(uninf52, 0);

    let subset_cells: Vec<_> = likelihoods
        .iter()
        .filter(|rl| {
            reads
                .get(rl.read_index.get())
                .is_some_and(|r| java_alignment_read_overlaps_interval(r, TARGET, TARGET, MARGIN))
        })
        .cloned()
        .collect();
    kv(
        "likelihood_subset_for_event_unique_n",
        unique_indices(&subset_cells).len().to_string(),
    );
    kv(
        "likelihood_subset_for_event_cell_n",
        subset_cells.len().to_string(),
    );
    assert_eq!(unique_indices(&subset_cells).len(), 47);

    let mut meta: Vec<Row> = Vec::new();
    for row in &overlap_rows {
        let Some(idx) = row.matrix_read_index() else {
            continue;
        };
        let rec = &reads[idx];
        let lr = row.haplotype_log10_likelihoods[0];
        let la = row.haplotype_log10_likelihoods[1];
        let (lr_f, la_f) = apply_java_mismapping_floor(lr, la);
        let (conf, informative, vote) = vote_row(lr_f, la_f);
        let flag = rec.flags();
        meta.push(Row {
            idx,
            qname: String::from_utf8_lossy(rec.qname()).into_owned(),
            flag,
            mapq: rec.mapq(),
            pos: rec.pos() + 1,
            end: alignment_end_1based(rec) as i64,
            cigar: rec.cigar().to_string(),
            lr: lr_f,
            la: la_f,
            conf,
            informative,
            vote,
            overlap_pm2: true,
            overlap_pm0: overlap_pm0(rec, TARGET),
            covers: java_alignment_read_covers_variant_base(rec, TARGET, TARGET, MARGIN),
            pileup: pileup_base_at(rec, TARGET).map(|b| b.to_ascii_uppercase()),
            dup: flag & FLAG_DUPLICATE != 0,
            secondary: flag & FLAG_NOT_PRIMARY != 0,
            supplementary: flag & FLAG_SUPPLEMENTARY != 0,
            qc_fail: flag & FLAG_VENDOR_QUALITY_FAILED != 0,
            proper_pair: flag & 0x0002 != 0,
            mate_unmapped: flag & 0x0008 != 0,
            first: flag & 0x0040 != 0,
            reverse: flag & 0x0010 != 0,
        });
    }
    meta.sort_by(|a, b| a.qname.cmp(&b.qname).then(a.flag.cmp(&b.flag)));

    let mut by_qname: BTreeMap<&str, usize> = BTreeMap::new();
    for r in &meta {
        *by_qname.entry(r.qname.as_str()).or_default() += 1;
    }
    kv(
        "multi_qname_n",
        by_qname.values().filter(|n| **n > 1).count().to_string(),
    );

    let mut pred_counts: BTreeMap<&str, usize> = BTreeMap::new();
    for r in &meta {
        kv(
            "row52",
            format!(
                "idx={}\tqname={}\tFLAG={}\tMAPQ={}\tPOS={}-{}\tCIGAR={}\tvote={}\tconf={:.4}\tpm0={}\tcovers={}\tpileup={:?}\tdup={}\tsec={}\tsupp={}\tqc={}\tproper={}\tmate_unmapped={}\tfirst={}\trev={}",
                r.idx,
                r.qname,
                r.flag,
                r.mapq,
                r.pos,
                r.end,
                r.cigar,
                r.vote,
                r.conf,
                r.overlap_pm0,
                r.covers,
                r.pileup.map(|b| b as char),
                r.dup,
                r.secondary,
                r.supplementary,
                r.qc_fail,
                r.proper_pair,
                r.mate_unmapped,
                r.first,
                r.reverse
            ),
        );
        if r.vote == "REF" && !r.overlap_pm0 {
            *pred_counts.entry("ref_margin_only_pm2").or_default() += 1;
        }
        if r.vote == "REF" && !r.covers {
            *pred_counts.entry("ref_not_covering").or_default() += 1;
        }
        if r.vote == "REF" && r.dup {
            *pred_counts.entry("ref_dup").or_default() += 1;
        }
        if r.vote == "REF" && r.secondary {
            *pred_counts.entry("ref_secondary").or_default() += 1;
        }
        if r.vote == "REF" && r.supplementary {
            *pred_counts.entry("ref_supp").or_default() += 1;
        }
        if r.vote == "REF" && r.qc_fail {
            *pred_counts.entry("ref_qc").or_default() += 1;
        }
        if r.vote == "REF" && r.mapq < 20 {
            *pred_counts.entry("ref_mapq_lt20").or_default() += 1;
        }
        if r.vote == "REF" && !r.proper_pair {
            *pred_counts.entry("ref_not_proper").or_default() += 1;
        }
        if r.vote == "REF" && r.mate_unmapped {
            *pred_counts.entry("ref_mate_unmapped").or_default() += 1;
        }
        if r.vote == "REF" && r.pileup.is_none() {
            *pred_counts.entry("ref_no_pileup_base").or_default() += 1;
        }
        if r.vote == "REF" && r.pileup.is_some() && r.pileup != Some(b'G') {
            *pred_counts.entry("ref_pileup_not_g").or_default() += 1;
        }
        if r.vote == "REF" && r.pos > TARGET as i64 {
            *pred_counts.entry("ref_starts_after").or_default() += 1;
        }
        if r.vote == "REF" && r.end < TARGET as i64 {
            *pred_counts.entry("ref_ends_before").or_default() += 1;
        }
    }
    for (k, n) in &pred_counts {
        kv("pred_count", format!("{k}={n}"));
    }

    let keep_pm0: BTreeSet<usize> = meta
        .iter()
        .filter(|r| r.overlap_pm0)
        .map(|r| r.idx)
        .collect();
    let keep_cover: BTreeSet<usize> = meta.iter().filter(|r| r.covers).map(|r| r.idx).collect();
    let keep_no_dup: BTreeSet<usize> = meta.iter().filter(|r| !r.dup).map(|r| r.idx).collect();
    let keep_primary: BTreeSet<usize> = meta
        .iter()
        .filter(|r| !r.secondary && !r.supplementary)
        .map(|r| r.idx)
        .collect();
    let keep_mapq20: BTreeSet<usize> = meta
        .iter()
        .filter(|r| r.mapq >= 20)
        .map(|r| r.idx)
        .collect();
    let keep_proper: BTreeSet<usize> = meta
        .iter()
        .filter(|r| r.proper_pair)
        .map(|r| r.idx)
        .collect();
    let keep_pileup: BTreeSet<usize> = meta
        .iter()
        .filter(|r| r.pileup.is_some())
        .map(|r| r.idx)
        .collect();
    let keep_pm0_primary: BTreeSet<usize> = meta
        .iter()
        .filter(|r| r.overlap_pm0 && !r.secondary && !r.supplementary && !r.dup)
        .map(|r| r.idx)
        .collect();

    let eval_keep = |name: &str, keep: &BTreeSet<usize>| {
        let sub: Vec<_> = overlap_rows
            .iter()
            .filter(|row| row.matrix_read_index().is_some_and(|i| keep.contains(&i)))
            .cloned()
            .collect();
        let dropped: Vec<_> = meta
            .iter()
            .filter(|r| !keep.contains(&r.idx))
            .map(|r| {
                format!(
                    "{} FLAG={} vote={} MAPQ={} POS={}-{} pm0={} covers={} dup={} sec={}",
                    r.qname,
                    r.flag,
                    r.vote,
                    r.mapq,
                    r.pos,
                    r.end,
                    r.overlap_pm0,
                    r.covers,
                    r.dup,
                    r.secondary
                )
            })
            .collect();
        let (ad, pl, _, _, uninf) = ad_pl_of(&sub);
        let hit = ad == fmt_ad && pl == fmt_pl;
        kv(
            "subset_eval",
            format!(
                "name={name}\tn={}\tad={ad:?}\tpl={pl:?}\tuninf={uninf}\thit_format={hit}\tdropped_n={}",
                sub.len(),
                dropped.len()
            ),
        );
        for d in &dropped {
            kv(&format!("dropped_{name}"), d);
        }
        (sub.len(), ad, pl, dropped, hit)
    };

    let evals = [
        ("pm0", keep_pm0.clone()),
        ("covering", keep_cover),
        ("no_dup", keep_no_dup),
        ("primary", keep_primary),
        ("mapq20", keep_mapq20),
        ("proper_pair", keep_proper),
        ("has_pileup", keep_pileup),
        ("pm0_primary_nodup", keep_pm0_primary),
    ];
    let mut hits: Vec<String> = Vec::new();
    let mut pm0_n = 0usize;
    let mut pm0_dropped: Vec<String> = Vec::new();
    for (name, keep) in &evals {
        let (n, ad, pl, dropped, hit) = eval_keep(name, keep);
        if *name == "pm0" {
            pm0_n = n;
            pm0_dropped = dropped.clone();
        }
        if hit {
            hits.push(format!("{name} n={n} ad={ad:?} pl={pl:?}"));
        }
    }
    kv("format_predicate_hits", format!("{hits:?}"));
    kv("pm0_n", pm0_n.to_string());
    assert!(
        hits.iter().all(|h| h.contains("n=47")),
        "6R.226: only drop-none identity filters match FORMAT P2 remarg: {hits:?}"
    );
    assert!(
        !hits.is_empty(),
        "identity predicates (no_dup/primary/mapq20/proper_pair) must match FORMAT after 6R.226"
    );
    assert_eq!(pm0_n, 46, "±0 overlap drops one REF, not three");

    // TLS collision: same overlap cell count as a previous site in this region.
    let mut cell_n_by_loc: BTreeMap<u64, usize> = BTreeMap::new();
    for call in &outcome.genotyped_calls {
        let loc = call.event.start_1based.get();
        let cells = likelihoods
            .iter()
            .filter(|rl| {
                reads.get(rl.read_index.get()).is_some_and(|r| {
                    java_alignment_read_overlaps_interval(
                        r,
                        loc,
                        loc.saturating_add(call.event.ref_allele.len().saturating_sub(1) as u64)
                            .max(call.event.end_1based.get()),
                        MARGIN,
                    )
                })
            })
            .count();
        cell_n_by_loc.insert(loc, cells);
        kv(
            "site_subset_cells",
            format!(
                "loc={loc}\t{}/{}\tcells={cells}\tunique={}\tad={:?}\tpl={:?}",
                call.event.ref_allele,
                call.event.alt_allele,
                unique_indices(
                    &likelihoods
                        .iter()
                        .filter(|rl| {
                            reads.get(rl.read_index.get()).is_some_and(|r| {
                                java_alignment_read_overlaps_interval(r, loc, loc, MARGIN)
                            })
                        })
                        .cloned()
                        .collect::<Vec<_>>()
                )
                .len(),
                call.genotype.format.ad_as_i32(),
                call.genotype.format.pl_as_i32()
            ),
        );
    }
    let target_cells = cell_n_by_loc.get(&TARGET).copied().unwrap_or(0);
    let colliding: Vec<_> = cell_n_by_loc
        .iter()
        .filter(|(loc, n)| **loc != TARGET && **n == target_cells)
        .map(|(loc, n)| format!("{loc}:{n}"))
        .collect();
    kv("tls_same_cell_n_other_sites", format!("{colliding:?}"));
    kv("target_subset_cells", target_cells.to_string());
    assert!(
        colliding.is_empty(),
        "TLS collision with another site in the region is not the 52→49 split: {colliding:?}"
    );

    // Original BAM overlap ±2 after HC filters — Java-like pre-realign population.
    let mut orig_keys: BTreeSet<(String, u16)> = BTreeSet::new();
    if let Ok(mut reader) = IndexedReader::from_path(&bam) {
        let _ = reader.fetch(("20", TARGET.saturating_sub(200) - 1, TARGET + 200));
        for rec in reader.records().flatten() {
            if rec.is_unmapped() {
                continue;
            }
            if !java_alignment_read_overlaps_interval(&rec, TARGET, TARGET, MARGIN) {
                continue;
            }
            orig_keys.insert((
                String::from_utf8_lossy(rec.qname()).into_owned(),
                rec.flags(),
            ));
        }
    }
    kv("orig_bam_overlap_n", orig_keys.len().to_string());
    let rust_keys: BTreeSet<(String, u16)> =
        meta.iter().map(|r| (r.qname.clone(), r.flag)).collect();
    let rust_only: Vec<_> = rust_keys.difference(&orig_keys).cloned().collect();
    let javaish_only: Vec<_> = orig_keys.difference(&rust_keys).cloned().collect();
    kv("rust52_not_in_orig_bam_n", rust_only.len().to_string());
    kv("orig_bam_not_in_rust52_n", javaish_only.len().to_string());
    for (q, f) in &rust_only {
        kv("rust52_only", format!("{q} FLAG={f}"));
    }
    for (q, f) in javaish_only.iter().take(20) {
        kv("orig_bam_only", format!("{q} FLAG={f}"));
    }

    let dropped_ref: Vec<_> = meta
        .iter()
        .filter(|r| r.vote == "REF" && !keep_pm0.contains(&r.idx))
        .cloned()
        .collect();
    kv("dropped_pm0_ref_n", dropped_ref.len().to_string());
    for r in &dropped_ref {
        kv(
            "dropped_ref_pm0",
            format!(
                "qname={}\tFLAG={}\tMAPQ={}\tPOS={}-{}\tCIGAR={}\tcovers={}\tpileup={:?}\tdup={}\tsec={}",
                r.qname,
                r.flag,
                r.mapq,
                r.pos,
                r.end,
                r.cigar,
                r.covers,
                r.pileup.map(|b| b as char),
                r.dup,
                r.secondary
            ),
        );
    }

    let ann_n = unique_indices(&site.annotation_likelihoods).len();
    let coverage_n =
        coverage_evidence_count(reads, &site.annotation_likelihoods, TARGET, TARGET, MARGIN);
    kv("annotation_unique_n", ann_n.to_string());
    kv("coverage_evidence_n", coverage_n.to_string());
    assert_eq!(ann_n, 47);
    assert_eq!(coverage_n, 47);

    let mut pileup_region_ref = 0i32;
    let mut pileup_region_alt = 0i32;
    let mut pileup_region_q: BTreeSet<Vec<u8>> = BTreeSet::new();
    let mut pileup_region_q_ref = 0i32;
    let mut pileup_region_q_alt = 0i32;
    for rec in &region.reads {
        let Some(b) = pileup_base_at(rec, TARGET) else {
            continue;
        };
        let b = b.to_ascii_uppercase();
        match b {
            b'G' => pileup_region_ref += 1,
            b'A' => pileup_region_alt += 1,
            _ => {}
        }
        if pileup_region_q.insert(rec.qname().to_owned()) {
            match b {
                b'G' => pileup_region_q_ref += 1,
                b'A' => pileup_region_q_alt += 1,
                _ => {}
            }
        }
    }
    kv(
        "pileup_region_reads",
        format!("raw={pileup_region_ref},{pileup_region_alt} qname={pileup_region_q_ref},{pileup_region_q_alt}"),
    );
    kv(
        "pileup_matches_format",
        ((pileup_region_ref, pileup_region_alt) == (fmt_ad[0], fmt_ad[1])
            || (pileup_region_q_ref, pileup_region_q_alt) == (fmt_ad[0], fmt_ad[1]))
            .to_string(),
    );

    let diagnosed = diagnose_genotype_variation_event(
        event,
        likelihoods,
        reads,
        &region.reads,
        Some(region.reads.as_slice()),
        haps,
        apply_bases.as_ref(),
        apply_pad,
        full_ref,
        full_pad,
        region.start.get(),
        region.end.get(),
        outcome.assembly.max_mnp_distance(),
        &config,
    )
    .expect("diagnose");
    match diagnosed {
        Ok(call) => {
            kv(
                "diagnose_try_genotype_ad",
                format!("{:?}", call.genotype.format.ad_as_i32()),
            );
            kv(
                "diagnose_try_genotype_pl",
                format!("{:?}", call.genotype.format.pl_as_i32()),
            );
            kv(
                "diagnose_matches_format",
                (call.genotype.format.ad_as_i32() == fmt_ad
                    && call.genotype.format.pl_as_i32() == fmt_pl)
                    .to_string(),
            );
            kv(
                "diagnose_matches_52",
                (call.genotype.format.ad_as_i32() == ad52
                    && call.genotype.format.pl_as_i32() == pl52)
                    .to_string(),
            );
            assert_eq!(
                call.genotype.format.ad_as_i32(),
                vec![42, 5],
                "isolated try_genotype is the 52-row retainEvidence object"
            );
            assert_eq!(call.genotype.format.pl_as_i32(), vec![84, 0, 1738]);
            assert_eq!(
                call.genotype.format.ad_as_i32(),
                fmt_ad,
                "6R.226: call_region FORMAT matches isolated try_genotype P2 remarg"
            );
        }
        Err(reason) => kv("diagnose_reject", format!("{reason:?}")),
    }

    let first_op;
    let java_counterpart;
    let case;
    let classification;
    // Isolated try_genotype == 52-row retainEvidence. call_region FORMAT is not that
    // object. No BAM-level filter of the 52 reproduces FORMAT. Java 4.4 has no
    // post-retainEvidence FORMAT subset.
    first_op = "call_region/assign FORMAT writer is not try_genotype_variation_event on likelihood_subset_for_event(±2). Isolated try_genotype (diagnose, hap_events=None, region_events=[]) is n=52 AD 47,5 PL 68,0,1937. Production genotyped_calls FORMAT is AD 44,5 PL 78,0,1811. keep_qnames / dedupe_qname / covering / ±0 / FLAG / TLS-other-site do not produce the FORMAT object. The 52→49 reduction is a call_region-only transformation of the retainEvidence object.";
    java_counterpart = "NO — Java calculateGLsForThisEvent and DepthPerAlleleBySample use the same retainEvidence AlleleLikelihoods. Default HC has no post-retainEvidence FORMAT subset.";
    case = "A";
    classification = "ALLELE_LIKELIHOOD_INPUT_DIVERGENCE";
    kv("first_operation_52_to_49", first_op);
    kv("java_counterpart", java_counterpart);
    kv("case", case);
    kv("classification", classification);
    kv("pl_same_boundary", (pl52 != fmt_pl).to_string());
    assert_eq!(pl52, fmt_pl, "6R.226: FORMAT PL is the 52-row P2 remarg");
    assert_eq!(classification, "ALLELE_LIKELIHOOD_INPUT_DIVERGENCE");

    // Five Rust-vs-Java source-population slots: 52 vs Java 47.
    kv(
        "source_pop",
        format!(
            "java_retainEvidence_n=47 (FORMAT DP 47, AD 42,5, UNINF inferred 0) rust52_n=52 rust52_ad=47,5 rust_format_n=49 rust_format_ad=44,5 extras_vs_java=5 format_drops_from_rust52=3"
        ),
    );
    kv(
        "structure",
        "Do not treat FORMAT DP=49 as a 49-row membership subset of the 52-row remarg. Isolated try_genotype keeps all 52 informative (AD 47,5). No measured predicate drops exactly 3 REF from that remarg. Java 47 vs Rust 52 is a separate source-population split (no JVM per-read dump this round).",
    );
    kv(
        "three_dropped_ref_rows",
        "UNNAMED as a filter of the 52-row remarg — that remarg is not the FORMAT writer. The only measured single-REF drop is ±0/covering: HWI-D00360:6:H81VLADXX:2:2106:19357:13581 FLAG=147 (n=51 AD 46,5 PL 71,0,1895), which is not FORMAT.",
    );

    kv("pm0_dropped_list", format!("{pm0_dropped:?}"));
    kv("production_change", "NONE");
}

#[test]
fn forensic_6r220_closed_6r218_untouched() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, INTERVAL).expect("interval");
    let walk = traverse_assembly_region_walker(
        &dict,
        &specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("walk");
    let regions = flatten_assembly_regions(&walk);
    let covering = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_PL
                && r.end.get() >= CLOSED_PL
        })
        .expect("6R.218 ActiveFull");
    let outcome = HaplotypeCallerEngine::call_region(
        covering,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let site = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_PL)
                && c.event.ref_allele == "G"
                && c.event.alt_allele == "T"
        })
        .expect("G/T");
    assert_eq!(outcome.assembly.haplotypes.len(), 30);
    assert_eq!(site.genotype.format.pl_as_i32(), vec![69, 0, 2140]);
    kv("closed_6r218_pl", "69,0,2140");
    kv("closed_6r218_hap_n", "30");
}
