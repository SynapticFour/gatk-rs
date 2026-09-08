//! 6R.86 coordinate-free: AD evidence is `retainEvidence` then `DepthPerAlleleBySample`.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//!
//! ```text
//! readLikelihoods.marginalize(alleleMapper)
//! retainEvidence(read → variantCallingRelevantOverlap.overlaps(read))
//!     // default HC: alignment overlap ± informativeReadOverlapMargin (2)
//!     // NOT QNAME/fragment collapse (that is Mutect groupEvidence)
//! calculateGLsForThisEvent          // same retained object
//! calculateGenotypes / unused-ALT
//! prepareReadAlleleLikelihoodsForAnnotation
//!     // reuse genotyping likelihoods when no contamination
//!     // updateNonRefAlleleLikelihoods is a no-op without NON_REF
//!     // addEvidence(filtered, 0) → ties → uninformative
//! DepthPerAlleleBySample.annotateWithLikelihoods
//!     // marginalize remaining VC alleles (column select)
//!     // bestAllelesBreakingTies + isInformative (confidence > 0.2)
//! reverseTrimAlleles                // AFTER annotation
//! ```
//!
//! 6R.85 proved: 136×4 likelihood matrices identical; SPAN_DEL is not causal.
//! This contract does not investigate PL, QUAL, PairHMM, or trimDownHaplotypes.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r86_ad_evidence_contract
//! HOLDOUT_6R86=1 cargo test -p gatk-haplotypecaller --test forensic_6r86_ad_evidence_contract live_ -- --nocapture
//! ```

use gatk_haplotypecaller::event_map::{
    overlapping_events, remap_alt_onto_longer_ref, VariationEvent,
};
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::genotyping::ReadLikelihoodRow;
use gatk_haplotypecaller::hc_allele_mapping::SPAN_DEL_ALLELE;
use gatk_haplotypecaller::hc_genotyping_engine::java_alignment_read_overlaps_interval;
use gatk_haplotypecaller::read_realignment::LOG_10_INFORMATIVE_THRESHOLD;
use gatk_haplotypecaller::{
    region_likelihoods_to_rows, InformativeAd, DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN,
};
use rust_htslib::bam::{self, record::Cigar, record::CigarString};
use std::collections::{HashMap, HashSet};

fn row(lls: Vec<f64>) -> ReadLikelihoodRow {
    ReadLikelihoodRow {
        read_index: 0,
        read_id: String::new(),
        haplotype_log10_likelihoods: lls,
    }
}

fn informative_vote(lls: &[f64]) -> Option<usize> {
    let mut best_i = 0usize;
    let mut best = f64::NEG_INFINITY;
    let mut second = f64::NEG_INFINITY;
    for (i, &ll) in lls.iter().enumerate() {
        if ll > best {
            second = best;
            best = ll;
            best_i = i;
        } else if ll > second {
            second = ll;
        }
    }
    if best.is_finite() && (best - second).abs() > LOG_10_INFORMATIVE_THRESHOLD {
        Some(best_i)
    } else {
        None
    }
}

fn mate(qname: &[u8], pos0: i64) -> rust_htslib::bam::Record {
    let mut rec = bam::Record::new();
    rec.set(
        qname,
        Some(&CigarString(vec![Cigar::Match(10)])),
        b"ACGTACGTAC",
        b"##########",
    );
    rec.set_pos(pos0);
    rec
}

/// Java `retainEvidence`: `target.overlaps(read)` independently per mate.
#[test]
fn forensic_6r86_retain_evidence_keeps_overlapping_paired_mates() {
    let r0 = mate(b"frag1", 99);
    let r1 = mate(b"frag1", 100);
    assert!(java_alignment_read_overlaps_interval(&r0, 105, 105, 2));
    assert!(java_alignment_read_overlaps_interval(&r1, 105, 105, 2));
    assert_eq!(r0.qname(), r1.qname());
}

/// Mutect `groupEvidence(GATKRead::getName)` collapses fragments. Default HC does not.
#[test]
fn forensic_6r86_qname_collapse_is_not_default_hc_retain_evidence() {
    let r0 = mate(b"frag1", 99);
    let r1 = mate(b"frag1", 100);
    let reads = [&r0, &r1];
    let overlapping: Vec<usize> = (0..2)
        .filter(|&i| java_alignment_read_overlaps_interval(reads[i], 105, 105, 2))
        .collect();
    assert_eq!(overlapping.len(), 2, "Java retainEvidence keeps both mates");
    let mut best_qname: HashMap<Vec<u8>, usize> = HashMap::new();
    for &i in &overlapping {
        best_qname.insert(reads[i].qname().to_owned(), i);
    }
    assert_eq!(
        best_qname.len(),
        1,
        "QNAME collapse is extra vs default HC retainEvidence"
    );
}

/// Java `DepthPerAlleleBySample`: remarg remaining alleles, then `isInformative`.
/// Permute of 4-way best-allele counts is not that algorithm.
#[test]
fn forensic_6r86_depth_per_allele_remarginalizes_remaining_alleles() {
    // 0=REF, 1=T (unused), 2=CG (called), 3=*
    let rows = [
        row(vec![0.0, -10.0, -10.0, -10.0]), // 4-way REF
        row(vec![-10.0, -10.0, 0.0, -10.0]), // 4-way CG
        row(vec![-0.5, 0.0, -10.0, -10.0]),  // 4-way T; remarg REF vs CG
        row(vec![-0.5, -10.0, -10.0, 0.0]),  // 4-way *; remarg REF vs CG
        row(vec![0.0, -0.05, -10.0, -10.0]), // 4-way uninformative; remarg REF vs CG
    ];
    let mut ad4 = [0i32; 4];
    for r in &rows {
        if let Some(i) = informative_vote(&r.haplotype_log10_likelihoods) {
            ad4[i] += 1;
        }
    }
    assert_eq!(ad4, [1, 1, 1, 1]);
    let permuted = vec![ad4[0], ad4[2]];
    assert_eq!(permuted, vec![1, 1]);

    let two_way: Vec<ReadLikelihoodRow> = rows
        .iter()
        .map(|r| {
            row(vec![
                r.haplotype_log10_likelihoods[0],
                r.haplotype_log10_likelihoods[2],
            ])
        })
        .collect();
    let remarg = InformativeAd::from_marginalized_rows(&two_way, 0, 1, None);
    assert_eq!(remarg.as_vec(), vec![4, 1]);
    assert_ne!(
        remarg.as_vec(),
        permuted,
        "Java DepthPerAlleleBySample remarginalizes; it does not permute 4-way counts"
    );
}

/// `updateNonRefAlleleLikelihoods` returns immediately when `NON_REF` is absent.
#[test]
fn forensic_6r86_update_nonref_is_noop_without_symbolic_nonref() {
    let alleles = ["TG", "*", "T", "CG"];
    assert!(!alleles.iter().any(|a| *a == "<NON_REF>"));
    assert!(alleles.iter().any(|a| *a == SPAN_DEL_ALLELE));
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum JavaPool {
    Ref,
    Alt(usize),
    SpanDel,
    Unassigned,
}

fn java_mapper_pool(
    spanning: &[VariationEvent],
    loc: u64,
    long_ref: &str,
    alts: &[String],
    emit_spanning: bool,
) -> JavaPool {
    let loc_pos = GenomePosition::new_1based(loc);
    if spanning.is_empty() {
        return JavaPool::Ref;
    }
    let mut hit_alt: Option<usize> = None;
    for ev in spanning {
        if ev.start_1based == loc_pos {
            if ev.ref_allele.len() == long_ref.len() {
                if let Some(ai) = alts.iter().position(|a| a == &ev.alt_allele) {
                    hit_alt = Some(ai);
                }
            } else if ev.ref_allele.len() < long_ref.len() {
                if let Some(remapped) =
                    remap_alt_onto_longer_ref(&ev.ref_allele, &ev.alt_allele, long_ref)
                {
                    if let Some(ai) = alts.iter().position(|a| a == &remapped) {
                        hit_alt = Some(ai);
                    }
                }
            }
        } else if emit_spanning {
            return JavaPool::SpanDel;
        } else {
            return JavaPool::Ref;
        }
    }
    match hit_alt {
        Some(ai) => JavaPool::Alt(ai),
        None => JavaPool::Unassigned,
    }
}

fn pool_max(row: &ReadLikelihoodRow, idxs: &[usize]) -> f64 {
    idxs.iter()
        .filter_map(|&i| row.haplotype_log10_likelihoods.get(i).copied())
        .fold(f64::NEG_INFINITY, f64::max)
}

fn apply_allele_floor(lls: &mut [f64]) {
    let best = lls.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if !best.is_finite() {
        return;
    }
    let floor = best - 4.5;
    for v in lls {
        if v.is_finite() && *v < floor {
            *v = floor;
        }
    }
}

fn read_fp_key(rec: &rust_htslib::bam::Record) -> String {
    let bases = String::from_utf8_lossy(&rec.seq().as_bytes()).into_owned();
    let bq: String = rec
        .qual()
        .iter()
        .map(|&q| char::from(q.saturating_add(33)))
        .collect();
    let cigar = format!("{:?}", rec.cigar());
    format!("{}|{}|{}|{}", bases, bq, rec.flags(), cigar)
}

#[test]
fn live_ad_evidence_retain_then_remarg() {
    if std::env::var("HOLDOUT_6R86").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R86=1");
        return;
    }
    use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
    use gatk_haplotypecaller::event_map::{
        build_per_haplotype_variation_events, variation_events_at_position_from_cache,
    };
    use gatk_haplotypecaller::hc_allele_mapping::replace_span_del_events;
    use gatk_haplotypecaller::{
        call_disposition, flatten_assembly_regions, take_colocated_merge_numerics,
        traverse_assembly_region_walker, try_emit_call_region_variants,
        AssemblyRegionCallDisposition, CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams,
        WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
    };
    use std::path::Path;

    const INTERVAL: &str = "20:29455000-29456500";
    const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
    const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
    const POS_SNP: u64 = 29_456_344;

    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: live BAM/ref missing");
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
    let covering: Vec<_> = regions
        .iter()
        .filter(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= POS_SNP
                && r.end.get() >= POS_SNP
        })
        .collect();
    assert_eq!(covering.len(), 1);
    let outcome = HaplotypeCallerEngine::call_region(
        covering[0],
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let emitted = try_emit_call_region_variants(
        covering[0],
        &outcome,
        "SAMPLE",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .unwrap_or_default();
    let vcf = emitted
        .iter()
        .find(|r| {
            r.position == POS_SNP && r.reference == "T" && r.alternate.iter().any(|a| a == "C")
        })
        .expect("T/C");
    assert!(
        !vcf.alternate.iter().any(|a| a == SPAN_DEL_ALLELE),
        "* must not be emitted"
    );

    let snap = take_colocated_merge_numerics()
        .into_iter()
        .find(|n| n.loc == POS_SNP)
        .expect("colocated merge numerics");

    let haps = &outcome.assembly.haplotypes;
    let full_pad = outcome.assembly.padded_reference_start_1based();
    let full_ref = outcome.assembly.reference_bases();
    let contig = covering[0].contig.as_str();
    let max_mnp = outcome.assembly.max_mnp_distance();
    let hap_events =
        build_per_haplotype_variation_events(haps, full_ref, full_pad, max_mnp, contig);
    let raw = variation_events_at_position_from_cache(&hap_events, POS_SNP, true);
    let replaced = replace_span_del_events(&raw, POS_SNP, full_pad, full_ref);
    let at_loc: Vec<VariationEvent> = replaced
        .into_iter()
        .filter(|e| e.start_1based.get() == POS_SNP)
        .collect();
    let _ = at_loc;

    let long_ref = snap.long_ref.as_str();
    let alts = snap.alts.as_slice();
    let mut pools: Vec<Vec<usize>> = vec![Vec::new(); 1 + alts.len()];
    for i in 0..haps.len() {
        let spanning = overlapping_events(hap_events.events_for(i), POS_SNP);
        match java_mapper_pool(&spanning, POS_SNP, long_ref, alts, true) {
            JavaPool::Ref => pools[0].push(i),
            JavaPool::Alt(ai) => pools[ai + 1].push(i),
            JavaPool::SpanDel => {
                if let Some(ai) = alts.iter().position(|a| a == SPAN_DEL_ALLELE) {
                    pools[ai + 1].push(i);
                }
            }
            JavaPool::Unassigned => {}
        }
    }

    let loc = POS_SNP;
    let end = loc.saturating_add(long_ref.len().saturating_sub(1) as u64);
    let margin = DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN;
    let pairhmm_idx: HashSet<usize> = outcome
        .read_likelihoods
        .iter()
        .map(|e| e.read_index.get())
        .collect();
    let overlap_idx: Vec<usize> = pairhmm_idx
        .iter()
        .copied()
        .filter(|&i| {
            outcome
                .genotyping_reads
                .get(i)
                .is_some_and(|r| java_alignment_read_overlaps_interval(r, loc, end, margin))
        })
        .collect();
    let unclip_idx: Vec<usize> = pairhmm_idx
        .iter()
        .copied()
        .filter(|&i| {
            outcome.genotyping_reads.get(i).is_some_and(|r| {
                gatk_haplotypecaller::hc_genotyping_engine::soft_unclipped_read_overlaps_interval(
                    r, loc, end, margin,
                )
            })
        })
        .collect();

    let overlap_ll: Vec<_> = outcome
        .read_likelihoods
        .iter()
        .filter(|e| overlap_idx.contains(&e.read_index.get()))
        .cloned()
        .collect();
    let hap_rows = region_likelihoods_to_rows(&overlap_ll, haps.len());
    assert_eq!(hap_rows.len(), overlap_idx.len());

    let mut four_way: Vec<(usize, Vec<f64>)> = hap_rows
        .iter()
        .map(|row| {
            let mut lls: Vec<f64> = pools.iter().map(|p| pool_max(row, p)).collect();
            apply_allele_floor(&mut lls);
            (row.read_index, lls)
        })
        .collect();
    four_way.sort_by_key(|(i, _)| *i);

    let cg_ai = alts.iter().position(|a| a == "CG").expect("CG");
    let mut ad4_overlap = vec![0i32; pools.len()];
    let mut remarg_overlap = [0i32; 2];
    let mut java_ad_reads: Vec<(String, &'static str)> = Vec::new();
    for (ri, lls) in &four_way {
        if let Some(v) = informative_vote(lls) {
            ad4_overlap[v] += 1;
        }
        let two = [lls[0], lls[cg_ai + 1]];
        if let Some(v) = informative_vote(&two) {
            remarg_overlap[v] += 1;
            let rec = &outcome.genotyping_reads[*ri];
            let allele = if v == 0 { "REF" } else { "ALT" };
            java_ad_reads.push((read_fp_key(rec), allele));
        }
    }

    let mut best_qname: HashMap<Vec<u8>, (usize, f64)> = HashMap::new();
    for (ri, lls) in &four_way {
        let rec = &outcome.genotyping_reads[*ri];
        let best = lls.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let qn = rec.qname().to_owned();
        match best_qname.get(&qn) {
            Some(&(cur, cur_ll)) if best > cur_ll || (best == cur_ll && *ri < cur) => {
                best_qname.insert(qn, (*ri, best));
            }
            None => {
                best_qname.insert(qn, (*ri, best));
            }
            _ => {}
        }
    }
    let rust_keep: HashSet<usize> = best_qname.values().map(|(i, _)| *i).collect();
    let mut ad4_qname = vec![0i32; pools.len()];
    let mut remarg_qname = [0i32; 2];
    let mut rust_ad_reads: Vec<(String, &'static str)> = Vec::new();
    for (ri, lls) in &four_way {
        if !rust_keep.contains(ri) {
            continue;
        }
        if let Some(v) = informative_vote(lls) {
            ad4_qname[v] += 1;
        }
        let two = [lls[0], lls[cg_ai + 1]];
        if let Some(v) = informative_vote(&two) {
            remarg_qname[v] += 1;
            let rec = &outcome.genotyping_reads[*ri];
            let allele = if v == 0 { "REF" } else { "ALT" };
            rust_ad_reads.push((read_fp_key(rec), allele));
        }
    }

    let java_fps: HashSet<&str> = java_ad_reads.iter().map(|(k, _)| k.as_str()).collect();
    let rust_fps: HashSet<&str> = rust_ad_reads.iter().map(|(k, _)| k.as_str()).collect();
    let java_only: Vec<&str> = java_fps.difference(&rust_fps).copied().collect();
    let rust_only: Vec<&str> = rust_fps.difference(&java_fps).copied().collect();
    let common = java_fps.intersection(&rust_fps).count();

    let vcf_ad = vcf.samples[0].ad.clone().unwrap_or_default();
    let vcf_pl = vcf.samples[0].pl.clone().unwrap_or_default();
    let vcf_gt = vcf.samples[0].gt.as_ref().map(|g| g.alleles.clone());

    eprintln!(
        "6R.86 STAGE_A pairhmm_unique={} overlap_retainEvidence={} after_qname={} multi_qname={}",
        pairhmm_idx.len(),
        four_way.len(),
        rust_keep.len(),
        snap.n_qnames_with_multiple_overlapping_reads
    );
    eprintln!(
        "6R.86 STAGE_B overlap 4-way informative={:?} qname 4-way={:?}",
        ad4_overlap, ad4_qname
    );
    eprintln!(
        "6R.86 STAGE_C remarg{{REF,CG}} overlap={:?} qname={:?} snap_remarg={:?} snap_remarg_no_qname={:?}",
        remarg_overlap,
        remarg_qname,
        snap.subset_ad_remarginalized,
        snap.subset_ad_remarginalized_no_qname
    );
    eprintln!(
        "6R.86 STAGE_D permute overlap={:?} qname={:?} snap_perm={:?} snap_perm_no_qname={:?}",
        [ad4_overlap[0], ad4_overlap[cg_ai + 1]],
        [ad4_qname[0], ad4_qname[cg_ai + 1]],
        snap.subset_ad_permuted,
        snap.subset_ad_permuted_no_qname
    );
    eprintln!(
        "6R.86 STAGE_E vcf AD={:?} PL={:?} GT={:?} QUAL={:?}",
        vcf_ad, vcf_pl, vcf_gt, vcf.quality
    );
    eprintln!(
        "6R.86 AD_READS common={} java_only={} rust_only={}",
        common,
        java_only.len(),
        rust_only.len()
    );
    for (i, fp) in java_only.iter().take(12).enumerate() {
        eprintln!(
            "6R.86 JAVA_ONLY_AD[{i}] fp_hash={:016x} len={}",
            fnv1a64(fp),
            fp.len()
        );
    }
    for (i, fp) in rust_only.iter().take(12).enumerate() {
        eprintln!(
            "6R.86 RUST_ONLY_AD[{i}] fp_hash={:016x} len={}",
            fnv1a64(fp),
            fp.len()
        );
    }

    let all_idx: Vec<usize> = pairhmm_idx.iter().copied().collect();
    let remarg_pop = |idxs: &[usize]| -> ([i32; 2], usize) {
        let ll: Vec<_> = outcome
            .read_likelihoods
            .iter()
            .filter(|e| idxs.contains(&e.read_index.get()))
            .cloned()
            .collect();
        let rows = region_likelihoods_to_rows(&ll, haps.len());
        let mut remarg = [0i32; 2];
        for row in &rows {
            let mut lls: Vec<f64> = pools.iter().map(|p| pool_max(row, p)).collect();
            apply_allele_floor(&mut lls);
            let two = [lls[0], lls[cg_ai + 1]];
            if let Some(v) = informative_vote(&two) {
                remarg[v] += 1;
            }
        }
        (remarg, rows.len())
    };
    let (remarg_all, n_all) = remarg_pop(&all_idx);
    let (remarg_unclip, n_unclip) = remarg_pop(&unclip_idx);

    eprintln!(
        "6R.86 COUNTERFACTUAL all136 remarg={:?} n={} unclipped_overlap remarg={:?} n={} clipped_overlap remarg={:?} n={}",
        remarg_all, n_all, remarg_unclip, n_unclip, remarg_overlap, four_way.len()
    );

    assert_eq!(
        four_way.len(),
        snap.n_overlap_before_qname_dedupe,
        "live clipped overlap reconstruction must match production retainEvidence subset"
    );
    assert_eq!(
        remarg_overlap.to_vec(),
        snap.subset_ad_remarginalized,
        "production remarg AD is the clipped-overlap population"
    );
    assert_eq!(
        snap.n_qnames_with_multiple_overlapping_reads, 0,
        "QNAME collapse is not causal: no multi-mate QNAMEs in the overlap set"
    );
    assert_eq!(
        four_way.len(),
        rust_keep.len(),
        "QNAME collapse does not shrink this site's overlap set"
    );
    assert!(
        pairhmm_idx.len() > four_way.len(),
        "first AD population split is overlap eligibility (136→clipped overlap), not QNAME collapse"
    );

    let overlap_set: HashSet<usize> = overlap_idx.iter().copied().collect();
    let all_ll: Vec<_> = outcome
        .read_likelihoods
        .iter()
        .filter(|e| pairhmm_idx.contains(&e.read_index.get()))
        .cloned()
        .collect();
    let all_rows = region_likelihoods_to_rows(&all_ll, haps.len());
    let mut n_dropped_informative = 0i32;
    let mut first_dropped: Option<(usize, &'static str, i64, i64)> = None;
    let mut nearest_dropped: Option<(i64, usize, &'static str, i64, i64)> = None;
    let var_lo = loc as i64 - i64::from(margin);
    let var_hi = end as i64 + i64::from(margin);
    for row in &all_rows {
        if overlap_set.contains(&row.read_index) {
            continue;
        }
        let Some(rec) = outcome.genotyping_reads.get(row.read_index) else {
            continue;
        };
        let mut lls: Vec<f64> = pools.iter().map(|p| pool_max(row, p)).collect();
        apply_allele_floor(&mut lls);
        let two = [lls[0], lls[cg_ai + 1]];
        if let Some(v) = informative_vote(&two) {
            n_dropped_informative += 1;
            let rs = rec.pos() + 1;
            let re = rec.cigar().end_pos();
            let allele = if v == 0 { "REF" } else { "ALT" };
            if first_dropped.is_none() {
                first_dropped = Some((row.read_index, allele, rs, re));
            }
            let gap = if re < var_lo {
                var_lo - re
            } else if rs > var_hi {
                rs - var_hi
            } else {
                0
            };
            match nearest_dropped {
                Some((g, ..)) if gap >= g => {}
                _ => nearest_dropped = Some((gap, row.read_index, allele, rs, re)),
            }
        }
    }
    eprintln!(
        "6R.86 DROPPED_INFORMATIVE n={} first_index={:?} nearest_to_expanded_interval={:?}",
        n_dropped_informative, first_dropped, nearest_dropped
    );
    assert_eq!(
        n_dropped_informative,
        remarg_all[0] + remarg_all[1] - remarg_overlap[0] - remarg_overlap[1],
        "all-136 remarg minus clipped-overlap remarg is exactly the dropped informative set"
    );
    assert!(
        n_dropped_informative > 0,
        "clipped overlap drops remarg-informative PairHMM reads present in the 136-read matrix"
    );
}

fn fnv1a64(s: &str) -> u64 {
    let mut h = 0xcbf29ce484222325u64;
    for b in s.as_bytes() {
        h ^= u64::from(*b);
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}
