//! 6R.219: first causal AD arrow at `20:29455379 G/A` (proof-only).
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE.
//!
//! Java FORMAT AD `42,5` vs Rust `44,5`. QUAL/GQ/PL/INFO are downstream.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r219_ad_membership_boundary -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{
    build_per_haplotype_variation_events, merged_alleles_for_genotyping,
    variation_events_at_position_from_cache,
};
use gatk_haplotypecaller::hc_allele_mapping::create_allele_mapper_with_events;
use gatk_haplotypecaller::hc_genotyping_engine::{
    alt_hap_indices_for_genotype_marginalization, biallelic_genotype_log10_likelihoods_gatk,
    java_alignment_read_covers_variant_base, java_alignment_read_overlaps_interval,
    marginalize_rows_to_biallelic_alleles, ref_hap_indices_for_genotype_marginalization,
    region_likelihoods_to_rows, HcGenotypingConfig, InformativeAd,
    DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN,
};
use gatk_haplotypecaller::query_index_at_reference_position;
use gatk_haplotypecaller::read_realignment::LOG_10_INFORMATIVE_THRESHOLD;
use gatk_haplotypecaller::read_unclip::alignment_end_1based;
use gatk_haplotypecaller::variant_site_hc_annotations::coverage_evidence_count;
use gatk_haplotypecaller::{
    begin_poorly_modeled_observe, call_disposition, flatten_assembly_regions,
    take_poorly_modeled_observe, traverse_assembly_region_walker, try_emit_call_region_variants,
    AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition, HaplotypeCallerEngine,
    ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use rust_htslib::bam::record::{Aux, CigarString};
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
const JAVA_INFO_DP: i32 = 47;
const JAVA_GQ: i32 = 84;
const JAVA_PL: [i32; 3] = [84, 0, 1738];
const JAVA_QUAL: f64 = 76.64;
const MARGIN: i32 = DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R219\t{key}\t{}", value.as_ref());
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

fn info_f64(info: &[InfoValue], key: &str) -> Option<f64> {
    for v in info {
        if let InfoValue::Float(k, xs) = v {
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
struct AdVote {
    idx: usize,
    qname: String,
    flag: u16,
    mapq: u8,
    pos: i64,
    end: i64,
    cigar: String,
    overlap: bool,
    lr: f64,
    la: f64,
    conf: f64,
    informative: bool,
    vote: &'static str,
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

fn orig_alignment(rec: &rust_htslib::bam::Record) -> (i64, i64, String, bool) {
    let has_orig = rec.aux(b"OP").is_ok() || rec.aux(b"OC").is_ok();
    let orig_pos0 = match rec.aux(b"OP") {
        Ok(Aux::I32(v)) => i64::from(v) - 1,
        Ok(Aux::U32(v)) => i64::from(v) - 1,
        Ok(Aux::I8(v)) => i64::from(v) - 1,
        Ok(Aux::U8(v)) => i64::from(v) - 1,
        Ok(Aux::I16(v)) => i64::from(v) - 1,
        Ok(Aux::U16(v)) => i64::from(v) - 1,
        _ => rec.pos(),
    };
    let orig_cigar = match rec.aux(b"OC") {
        Ok(Aux::String(s)) => s.to_string(),
        _ => rec.cigar().to_string(),
    };
    let ref_len = orig_cigar_ref_len(&orig_cigar).max(1);
    let start = orig_pos0 + 1;
    let end = orig_pos0 + ref_len;
    (start, end, orig_cigar, has_orig)
}

fn orig_cigar_ref_len(cigar: &str) -> i64 {
    let mut n = 0i64;
    let mut acc = 0i64;
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            acc = acc * 10 + i64::from(ch as u8 - b'0');
            continue;
        }
        if matches!(ch, 'M' | 'D' | 'N' | '=' | 'X') {
            n += acc;
        }
        acc = 0;
    }
    n
}

fn orig_overlaps(start: i64, end: i64, loc: u64, margin: i32) -> bool {
    let m = i64::from(margin.max(0));
    let loc = loc as i64;
    loc - m <= end && start <= loc + m
}

#[test]
fn forensic_6r219_java_ad_contract_pin() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE — proof-only");
    kv(
        "java_ad_lifecycle",
        "assignGenotypeLikelihoods: marginalize → retainEvidence(±2 overlap) → calculateGLsForThisEvent (PL only) → prepareReadAlleleLikelihoodsForAnnotation (reuse genotyping AlleleLikelihoods; addEvidence(filtered,0) ties) → DepthPerAlleleBySample.annotateWithLikelihoods: identity remarg remaining vc.getAlleles() → bestAllelesBreakingTies → isInformative(confidence > 0.2) → FORMAT AD",
    );
    kv(
        "java_ad_source",
        "DepthPerAlleleBySample on remaining call alleles of the retainEvidence AlleleLikelihoods (same object as calculateGLs when contamination off)",
    );
    assert_eq!(LOG_10_INFORMATIVE_THRESHOLD, 0.2);
    assert_eq!(MARGIN, 2);
}

#[test]
fn forensic_6r219_live_ad_membership_boundary() {
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

    begin_poorly_modeled_observe();
    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let poorly = take_poorly_modeled_observe();

    let haps = &outcome.assembly.haplotypes;
    let reads = &outcome.genotyping_reads;
    let likelihoods = &outcome.read_likelihoods;
    kv("hap_n", haps.len().to_string());
    kv("genotyping_read_n", reads.len().to_string());
    kv("pairhmm_cell_n", likelihoods.len().to_string());
    kv(
        "pairhmm_unique_reads",
        unique_indices(likelihoods).len().to_string(),
    );

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
    kv(
        "cache_events",
        cache_only
            .iter()
            .map(|e| format!("{}/{}", e.ref_allele, e.alt_allele))
            .collect::<Vec<_>>()
            .join(";"),
    );
    kv(
        "merged_alleles",
        format!("{:?}", merged_alleles_for_genotyping(&cache_only, TARGET)),
    );
    assert_eq!(cache_only.len(), 1, "biallelic EventMap at target");
    assert_eq!(cache_only[0].ref_allele, TARGET_REF);
    assert_eq!(cache_only[0].alt_allele, TARGET_ALT);

    let apply_bases = outcome.assembly.apply_bases_shared();
    let apply_pad = haps
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.map(|g| g.start_1based()))
        .unwrap_or(full_pad);
    let mapping = create_allele_mapper_with_events(
        &cache_only[0],
        TARGET,
        haps,
        apply_pad,
        apply_bases.as_ref(),
        outcome.assembly.max_mnp_distance(),
        emit_spanning,
        Some(&hap_events),
    );
    kv(
        "mapper",
        format!(
            "ref={}\talt={}",
            mapping.ref_haplotype_indices.len(),
            mapping.alt_haplotype_indices.len()
        ),
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
    kv(
        "extra_alt_alleles",
        if site.extra_alt_alleles.is_empty() {
            "(none)".into()
        } else {
            site.extra_alt_alleles.join(",")
        },
    );
    kv(
        "post_merge_unused_alt_subset",
        site.post_merge_unused_alt_subset.to_string(),
    );
    assert!(site.extra_alt_alleles.is_empty());
    assert!(!site.post_merge_unused_alt_subset);
    let pl = site.genotype.format.pl_as_i32();
    let ad = site.genotype.format.ad_as_i32();
    kv("rust_gt_pl", format!("{pl:?}"));
    kv("rust_gt_ad", format!("{ad:?}"));
    kv("rust_gt_gq", site.genotype.format.gq.as_i32().to_string());
    kv("rust_gt_dp", site.genotype.format.dp.as_i32().to_string());
    kv("java_ad", "42,5");
    kv("java_fmt_dp", JAVA_FMT_DP.to_string());
    kv("java_pl", "84,0,1738");
    assert_eq!(
        ad, JAVA_AD,
        "6R.239: FORMAT AD matches Java after sentinel fix"
    );

    let emitted =
        try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let rec = emitted
        .iter()
        .find(|r| r.position == TARGET && r.reference == TARGET_REF)
        .expect("emit G/A");
    let s = &rec.samples[0];
    kv(
        "vcf_gt",
        format!("{:?}", s.gt.as_ref().map(|g| g.to_string())),
    );
    kv("vcf_ad", format!("{:?}", s.ad.as_deref()));
    kv("vcf_dp", format!("{:?}", s.dp));
    kv("vcf_gq", format!("{:?}", s.gq));
    kv("vcf_pl", format!("{:?}", s.pl.as_deref()));
    kv("vcf_qual", format!("{:?}", rec.quality));
    kv("info_dp", format!("{:?}", info_i32(&rec.info, "DP")));
    kv("info_qd", format!("{:?}", info_f64(&rec.info, "QD")));
    kv("info_mq", format!("{:?}", info_f64(&rec.info, "MQ")));
    kv("info_fs", format!("{:?}", info_f64(&rec.info, "FS")));
    kv("info_sor", format!("{:?}", info_f64(&rec.info, "SOR")));
    assert_eq!(s.ad.as_deref(), Some([42, 5].as_slice()));
    assert_eq!(s.dp, Some(47));
    assert_eq!(info_i32(&rec.info, "DP"), Some(47));

    let config = HcGenotypingConfig::strict_java();
    let ref_hap = haps.iter().find(|h| h.is_reference).expect("ref hap");
    let ref_bytes = apply_bases.as_ref();
    let ref_pool =
        ref_hap_indices_for_genotype_marginalization(&mapping, haps, &config, Some(&cache_only[0]));
    let alt_pool = alt_hap_indices_for_genotype_marginalization(
        &mapping,
        haps,
        &cache_only[0],
        ref_hap,
        apply_pad,
        ref_bytes,
        outcome.assembly.max_mnp_distance(),
        &region.contig,
        &config,
    );
    kv(
        "pools",
        format!("ref={}\talt={}", ref_pool.len(), alt_pool.len()),
    );

    let rows = region_likelihoods_to_rows(likelihoods, haps.len());
    let marg_all = marginalize_rows_to_biallelic_alleles(&rows, &ref_pool, &alt_pool);
    let ad_all = InformativeAd::from_marginalized_rows(&marg_all, 0, 1, None);
    kv("ad_all_pairhmm_unique", format!("{ad_all:?}"));

    let overlap_idx: BTreeSet<usize> = reads
        .iter()
        .enumerate()
        .filter(|(_, rec)| java_alignment_read_overlaps_interval(rec, TARGET, TARGET, MARGIN))
        .map(|(i, _)| i)
        .collect();
    kv("overlap_genotyping_reads", overlap_idx.len().to_string());

    let overlap_rows: Vec<_> = marg_all
        .iter()
        .filter(|row| {
            row.matrix_read_index()
                .is_some_and(|i| overlap_idx.contains(&i))
        })
        .cloned()
        .collect();
    let ad_overlap = InformativeAd::from_marginalized_rows(&overlap_rows, 0, 1, None);
    kv("ad_retainEvidence_overlap", format!("{ad_overlap:?}"));
    kv(
        "n_retainEvidence_overlap_rows",
        overlap_rows.len().to_string(),
    );

    let ann_idx = unique_indices(&site.annotation_likelihoods);
    kv("annotation_unique_n", ann_idx.len().to_string());
    let coverage_n =
        coverage_evidence_count(reads, &site.annotation_likelihoods, TARGET, TARGET, MARGIN);
    kv("coverage_evidence_n", coverage_n.to_string());
    let ann_rows: Vec<_> = marg_all
        .iter()
        .filter(|row| {
            row.matrix_read_index()
                .is_some_and(|i| ann_idx.contains(&i))
        })
        .cloned()
        .collect();
    let ad_ann = InformativeAd::from_marginalized_rows(&ann_rows, 0, 1, None);
    kv("ad_annotation_object", format!("{ad_ann:?}"));

    let pairhmm_idx = unique_indices(likelihoods);
    kv(
        "overlap_not_in_pairhmm",
        overlap_idx.difference(&pairhmm_idx).count().to_string(),
    );
    kv(
        "pairhmm_overlap_not_in_ann",
        overlap_idx
            .intersection(&pairhmm_idx)
            .cloned()
            .collect::<BTreeSet<_>>()
            .difference(&ann_idx)
            .count()
            .to_string(),
    );

    let mut votes: Vec<AdVote> = Vec::new();
    let mut by_qname: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for row in &overlap_rows {
        let Some(idx) = row.matrix_read_index() else {
            continue;
        };
        let rec = &reads[idx];
        let lr = row.haplotype_log10_likelihoods[0];
        let la = row.haplotype_log10_likelihoods[1];
        let (conf, informative, vote) = vote_row(lr, la);
        let qname = String::from_utf8_lossy(rec.qname()).into_owned();
        by_qname.entry(qname.clone()).or_default().push(idx);
        let cigar = rec.cigar().to_string();
        votes.push(AdVote {
            idx,
            qname,
            flag: rec.flags(),
            mapq: rec.mapq(),
            pos: rec.pos() + 1,
            end: alignment_end_1based(rec) as i64,
            cigar,
            overlap: true,
            lr,
            la,
            conf,
            informative,
            vote,
        });
    }
    votes.sort_by(|a, b| a.qname.cmp(&b.qname).then(a.flag.cmp(&b.flag)));

    let mut ref_n = 0i32;
    let mut alt_n = 0i32;
    let mut uninf_n = 0i32;
    for v in &votes {
        match v.vote {
            "REF" => ref_n += 1,
            "ALT" => alt_n += 1,
            _ => uninf_n += 1,
        }
        kv(
            "ad_row",
            format!(
                "idx={}\tqname={}\tFLAG={}\tMAPQ={}\tPOS={}-{}\tCIGAR={}\tlr={:.6}\tla={:.6}\tconf={:.6}\tinf={}\tvote={}",
                v.idx,
                v.qname,
                v.flag,
                v.mapq,
                v.pos,
                v.end,
                v.cigar,
                v.lr,
                v.la,
                v.conf,
                v.informative,
                v.vote
            ),
        );
    }
    kv(
        "vote_counts",
        format!(
            "REF={ref_n}\tALT={alt_n}\tUNINF={uninf_n}\tn={}",
            votes.len()
        ),
    );
    kv(
        "multi_qname_n",
        by_qname
            .values()
            .filter(|xs| xs.len() > 1)
            .count()
            .to_string(),
    );
    kv(
        "format_vs_retainEvidence_remarg",
        format!("FORMAT={ad:?} retainEvidence_overlap={ad_overlap:?}"),
    );
    assert_eq!(
        ad_overlap.as_vec(),
        vec![ref_n, alt_n],
        "vote dump must match InformativeAd on the overlap remarg object"
    );
    assert_eq!(
        ad_overlap.as_vec(),
        ad,
        "6R.226: FORMAT AD is retainEvidence overlap informative remarg (P2)"
    );

    let subset_cells: Vec<_> = likelihoods
        .iter()
        .filter(|rl| {
            reads
                .get(rl.read_index.get())
                .is_some_and(|r| java_alignment_read_overlaps_interval(r, TARGET, TARGET, MARGIN))
        })
        .cloned()
        .collect();
    let subset_rows = region_likelihoods_to_rows(&subset_cells, haps.len());
    let subset_marg = marginalize_rows_to_biallelic_alleles(&subset_rows, &ref_pool, &alt_pool);
    let ad_subset = InformativeAd::from_marginalized_rows(&subset_marg, 0, 1, None);
    let pl_subset = gls_to_pl(&biallelic_genotype_log10_likelihoods_gatk(
        &subset_marg,
        0,
        1,
    ));
    kv(
        "site_score_subset_unique_n",
        unique_indices(&subset_cells).len().to_string(),
    );
    kv("ad_from_subset_cells", format!("{ad_subset:?}"));
    kv("pl_from_subset_cells", format!("{pl_subset:?}"));
    let gls_overlap = biallelic_genotype_log10_likelihoods_gatk(&overlap_rows, 0, 1);
    let pl_overlap = gls_to_pl(&gls_overlap);
    kv("pl_from_retainEvidence_overlap", format!("{pl_overlap:?}"));
    kv("pl_production", format!("{pl:?}"));

    let mut floored_rows = overlap_rows.clone();
    for row in &mut floored_rows {
        let (lr, la) = apply_java_mismapping_floor(
            row.haplotype_log10_likelihoods[0],
            row.haplotype_log10_likelihoods[1],
        );
        row.haplotype_log10_likelihoods[0] = lr;
        row.haplotype_log10_likelihoods[1] = la;
    }
    let ad_floor = InformativeAd::from_marginalized_rows(&floored_rows, 0, 1, None);
    let pl_floor = gls_to_pl(&biallelic_genotype_log10_likelihoods_gatk(
        &floored_rows,
        0,
        1,
    ));
    kv("ad_after_mismapping_floor", format!("{ad_floor:?}"));
    kv("pl_after_mismapping_floor", format!("{pl_floor:?}"));
    kv(
        "site_score_object_matches_overlap",
        (pl_floor == pl).to_string(),
    );
    if pl_floor != pl {
        let ref_idxs: Vec<usize> = floored_rows
            .iter()
            .enumerate()
            .filter(|(_, row)| {
                let lr = row.haplotype_log10_likelihoods[0];
                let la = row.haplotype_log10_likelihoods[1];
                vote_row(lr, la).2 == "REF"
            })
            .map(|(i, _)| i)
            .collect();
        let mut single_hits = Vec::new();
        for &drop in &ref_idxs {
            let sub: Vec<_> = floored_rows
                .iter()
                .enumerate()
                .filter(|(i, _)| *i != drop)
                .map(|(_, r)| r.clone())
                .collect();
            let p = gls_to_pl(&biallelic_genotype_log10_likelihoods_gatk(&sub, 0, 1));
            if p == pl {
                let idx = floored_rows[drop].matrix_read_index().unwrap_or(usize::MAX);
                single_hits.push(format!("drop_idx={idx}"));
            }
        }
        kv("pl_match_drop_one_ref", format!("{single_hits:?}"));
        let mut triple_hits = 0u32;
        let mut first_triple = String::new();
        for a in 0..ref_idxs.len() {
            for b in (a + 1)..ref_idxs.len() {
                for c in (b + 1)..ref_idxs.len() {
                    let da = ref_idxs[a];
                    let db = ref_idxs[b];
                    let dc = ref_idxs[c];
                    let sub: Vec<_> = floored_rows
                        .iter()
                        .enumerate()
                        .filter(|(i, _)| *i != da && *i != db && *i != dc)
                        .map(|(_, r)| r.clone())
                        .collect();
                    let p = gls_to_pl(&biallelic_genotype_log10_likelihoods_gatk(&sub, 0, 1));
                    if p == pl {
                        triple_hits += 1;
                        if first_triple.is_empty() {
                            let ia = floored_rows[da].matrix_read_index().unwrap_or(usize::MAX);
                            let ib = floored_rows[db].matrix_read_index().unwrap_or(usize::MAX);
                            let ic = floored_rows[dc].matrix_read_index().unwrap_or(usize::MAX);
                            first_triple = format!("{ia},{ib},{ic}");
                        }
                    }
                }
            }
        }
        kv("pl_match_drop_three_ref_n", triple_hits.to_string());
        kv("pl_match_drop_three_ref_first", first_triple);
    }

    let mut pileup_ref = 0i32;
    let mut pileup_alt = 0i32;
    let mut pileup_other = 0i32;
    let mut pileup_qnames: BTreeSet<Vec<u8>> = BTreeSet::new();
    let mut pileup_qname_ref = 0i32;
    let mut pileup_qname_alt = 0i32;
    for rec in reads {
        let Some(b) = pileup_base_at(rec, TARGET) else {
            continue;
        };
        let b = b.to_ascii_uppercase();
        match b {
            b'G' => pileup_ref += 1,
            b'A' => pileup_alt += 1,
            _ => pileup_other += 1,
        }
        if pileup_qnames.insert(rec.qname().to_owned()) {
            match b {
                b'G' => pileup_qname_ref += 1,
                b'A' => pileup_qname_alt += 1,
                _ => {}
            }
        }
    }
    let pileup_raw = (pileup_ref, pileup_alt);
    let pileup_qname = (pileup_qname_ref, pileup_qname_alt);
    kv(
        "pileup_ad_raw",
        format!("{pileup_raw:?} other={pileup_other}"),
    );
    kv("pileup_ad_qname_dedupe", format!("{pileup_qname:?}"));

    let cover_rows: Vec<_> = overlap_rows
        .iter()
        .filter(|row| {
            row.matrix_read_index().is_some_and(|i| {
                reads.get(i).is_some_and(|r| {
                    java_alignment_read_covers_variant_base(r, TARGET, TARGET, MARGIN)
                })
            })
        })
        .cloned()
        .collect();
    let ad_cover = InformativeAd::from_marginalized_rows(&cover_rows, 0, 1, None);
    kv("ad_covering_base_overlap", format!("{ad_cover:?}"));
    kv("n_covering_base_rows", cover_rows.len().to_string());

    let mut orig_overlap_rows = Vec::new();
    let mut orig_only = Vec::new();
    let mut realign_only = Vec::new();
    let mut pileup_skip_ref = Vec::new();
    for v in &votes {
        let rec = &reads[v.idx];
        let (os, oe, oc, has_orig) = orig_alignment(rec);
        let orig_ov = orig_overlaps(os, oe, TARGET, MARGIN);
        let covers = java_alignment_read_covers_variant_base(rec, TARGET, TARGET, MARGIN);
        let pb = pileup_base_at(rec, TARGET);
        kv(
            "ad_geom",
            format!(
                "idx={}\tqname={}\tvote={}\tcovers={}\tpileup_base={:?}\thas_OC/OP={}\torig={}-{} {}\torig_overlap={}",
                v.idx,
                v.qname,
                v.vote,
                covers,
                pb.map(|b| b as char),
                has_orig,
                os,
                oe,
                oc,
                orig_ov
            ),
        );
        if orig_ov {
            orig_overlap_rows.push(v.idx);
        }
        if orig_ov && !v.overlap {
            orig_only.push(v.qname.clone());
        }
        if !orig_ov && has_orig {
            realign_only.push(format!("{} FLAG={} vote={}", v.qname, v.flag, v.vote));
        }
        if v.vote == "REF" && (pb.is_none() || pb.map(|b| b.to_ascii_uppercase()) != Some(b'G')) {
            pileup_skip_ref.push(format!(
                "{} FLAG={} POS={}-{} CIGAR={} pileup={:?} covers={}",
                v.qname,
                v.flag,
                v.pos,
                v.end,
                v.cigar,
                pb.map(|b| b as char),
                covers
            ));
        }
    }
    kv("orig_only_n", orig_only.len().to_string());
    kv("realign_only_overlap_n", realign_only.len().to_string());
    for row in &realign_only {
        kv("realign_only_overlap", row);
    }
    kv(
        "remarg_ref_not_in_pileup_n",
        pileup_skip_ref.len().to_string(),
    );
    for row in &pileup_skip_ref {
        kv("remarg_ref_not_in_pileup", row);
    }

    let dropped_pm: Vec<_> = poorly
        .iter()
        .filter(|r| !r.rust_keep)
        .map(|r| {
            format!(
                "{} FLAG={} keep_java={} extra_retain={}",
                r.qname, r.flags, r.java_equiv_keep, r.extra_retain
            )
        })
        .collect();
    kv("poorly_modeled_dropped_n", dropped_pm.len().to_string());
    for row in dropped_pm.iter().take(12) {
        kv("poorly_modeled_drop", row);
    }
    let extra_retain: Vec<_> = poorly
        .iter()
        .filter(|r| r.extra_retain)
        .map(|r| format!("{} FLAG={}", r.qname, r.flags))
        .collect();
    kv(
        "poorly_modeled_extra_retain_n",
        extra_retain.len().to_string(),
    );

    let fmt_dp = site.genotype.format.dp.as_i32();
    kv(
        "dp_vs_ad",
        format!(
            "fmt_dp={fmt_dp} format_ad_sum={} remarg_ad_sum={} java_fmt_dp={JAVA_FMT_DP} java_ad_sum={} info_dp={:?} java_info_dp={JAVA_INFO_DP} pileup_raw={pileup_raw:?} pileup_qname={pileup_qname:?}",
            ad[0] + ad[1],
            ref_n + alt_n,
            JAVA_AD[0] + JAVA_AD[1],
            info_i32(&rec.info, "DP")
        ),
    );
    assert_eq!(fmt_dp, ad[0] + ad[1], "FORMAT DP tracks FORMAT AD sum");
    assert_eq!(fmt_dp, JAVA_FMT_DP, "6R.239: FORMAT DP matches Java");
    assert_eq!(
        info_i32(&rec.info, "DP"),
        Some(fmt_dp),
        "6R.226: INFO DP tracks the P2 AD-sum FORMAT DP on this site"
    );

    let first_arrow;
    let classification;
    if ad != JAVA_AD.to_vec() && ad_overlap.as_vec() != ad && ad_ann.as_vec() != ad {
        first_arrow = "FORMAT AD is a FORMAT-specific genotyping subset (44,5 / DP=49), not Java DepthPerAlleleBySample on retainEvidence AlleleLikelihoods (overlap remarg 47,5 / n=52) and not the annotation object";
        classification = "ALLELE_LIKELIHOOD_INPUT_DIVERGENCE";
    } else if ad_ann.as_vec() == JAVA_AD.to_vec() && ad != JAVA_AD.to_vec() {
        first_arrow = "FORMAT AD is not the annotation AlleleLikelihoods object; annotation remarg matches Java 42,5";
        classification = "ALLELE_LIKELIHOOD_INPUT_DIVERGENCE";
    } else if pileup_raw == (ad[0], ad[1]) && ad_overlap.as_vec() != ad {
        first_arrow = "FORMAT AD equals pileup covering-base counts, not DepthPerAlleleBySample informative remarg of retainEvidence overlap";
        classification = "AD_AGGREGATION_DIVERGENCE";
    } else if ad_cover.as_vec() == ad && ad_overlap.as_vec() != ad {
        first_arrow = "FORMAT AD equals covering-base overlap remarg; Java DepthPerAlleleBySample uses retainEvidence overlap, not covering-base";
        classification = "READ_MEMBERSHIP_DIVERGENCE";
    } else if ad_overlap.as_vec() != ad {
        first_arrow = "FORMAT AD is not retainEvidence overlap informative remarg";
        classification = "AD_AGGREGATION_DIVERGENCE";
    } else if uninf_n > 0 {
        first_arrow = "same overlap set; informative/best-allele filter admits extra REF votes";
        classification = "INFORMATIVE_EVIDENCE_DIVERGENCE";
    } else {
        first_arrow =
            "6R.239: retainEvidence overlap informative remarg FORMAT AD matches Java 42,5";
        classification = "NO_DIVERGENCE_AT_THIS_BOUNDARY";
    }
    kv("first_divergent_arrow", first_arrow);
    kv("classification", classification);
    assert_eq!(classification, "NO_DIVERGENCE_AT_THIS_BOUNDARY");
    kv(
        "pairhmm_causal",
        "winner identities dumped per overlap row; residuals not assumed causal",
    );
    kv("kbest_policy", "legacy_1024 unchanged");
    kv("production_change", "NONE");

    let extra_ref = ad[0] - JAVA_AD[0];
    kv("extra_ref_votes_format_vs_java", extra_ref.to_string());
    kv(
        "extra_ref_votes_remarg_vs_java",
        (ref_n - JAVA_AD[0]).to_string(),
    );
    kv(
        "extra_ref_votes_remarg_vs_format",
        (ref_n - ad[0]).to_string(),
    );
    assert_eq!(extra_ref, 0, "6R.239: P2 remarg matches Java REF votes");
    assert_eq!(ad[1], JAVA_AD[1], "ALT AD already matches Java");
    let _ = (JAVA_GQ, JAVA_PL, JAVA_QUAL);
}

#[test]
fn forensic_6r219_closed_6r218_untouched() {
    kv(
        "guard",
        "6R.218 cycle abort stays; this round does not reopen 20:29455015",
    );
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
    assert_eq!(site.genotype.format.gq.as_i32(), 69);
    kv("closed_6r218_pl", "69,0,2140");
    kv("closed_6r218_hap_n", "30");
}
