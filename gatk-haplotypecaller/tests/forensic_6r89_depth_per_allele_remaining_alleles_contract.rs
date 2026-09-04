//! 6R.89 coordinate-free: DepthPerAlleleBySample remarg + bestAllelesBreakingTies + isInformative.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//!
//! ```text
//! DepthPerAlleleBySample.annotate
//!   alleles = LinkedHashSet(vc.getAlleles())          // remaining call alleles, identity
//!   likelihoods.alleles().containsAll(alleles)
//! annotateWithLikelihoods
//!   alleleSubset = {allele → [allele]}                // column select, NOT log-sum
//!   subsetted = likelihoods.marginalize(alleleSubset) // max over mapped old columns
//!   bestAllelesBreakingTies(sample)                   // REF priority 1.0 vs alt 0
//!     searchBestAllele(..., canBeReference=true, priorities)
//!     tie-break iff (best - second) < 0.2
//!     among alleles within 0.2 of best, pick highest priority (REF)
//!     confidence = best==second ? 0 : best - second
//!   isInformative: confidence > 0.2                   // strict greater-than
//!   count by allele identity into vc allele order
//! ```
//!
//! 6R.88 proved the 62×4 AD input object is equivalent. This round starts at remarg.
//! PL, QUAL, PairHMM, overlap, and permute-vs-remarg production patches are out of scope.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r89_depth_per_allele_remaining_alleles_contract
//! HOLDOUT_6R89=1 cargo test -p gatk-haplotypecaller --test forensic_6r89_depth_per_allele_remaining_alleles_contract live_ -- --nocapture
//! ```

use gatk_haplotypecaller::event_map::{
    overlapping_events, remap_alt_onto_longer_ref, VariationEvent,
};
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::genotyping::ReadLikelihoodRow;
use gatk_haplotypecaller::hc_allele_mapping::SPAN_DEL_ALLELE;
use gatk_haplotypecaller::hc_genotyping_engine::java_alignment_read_overlaps_interval;
use gatk_haplotypecaller::read_realignment::LOG_10_INFORMATIVE_THRESHOLD;
use gatk_haplotypecaller::{region_likelihoods_to_rows, DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN};

const JAVA_INFORMATIVE: f64 = 0.2;

/// Java `AlleleLikelihoods.marginalize` with identity map: per new allele, max of mapped old columns.
fn remarg_identity(old: &[f64], keep: &[usize]) -> Vec<f64> {
    keep.iter()
        .map(|&i| old.get(i).copied().unwrap_or(f64::NEG_INFINITY))
        .collect()
}

/// Java `searchBestAllele` + default REF tie-break (`priority = isReference ? 1 : 0`).
/// `ref_i` is the reference column in `lls` (remaining-allele space).
fn java_best_breaking_ties(lls: &[f64], ref_i: usize) -> (usize, usize, f64, bool) {
    let n = lls.len();
    if n == 0 {
        return (0, 0, 0.0, false);
    }
    let mut best_i = 0usize;
    let mut second_i = 0usize;
    let mut best = lls[0];
    let mut second = f64::NEG_INFINITY;
    for a in 1..n {
        let cand = lls[a];
        if cand > best {
            second_i = best_i;
            second = best;
            best_i = a;
            best = cand;
        } else if cand > second {
            second_i = a;
            second = cand;
        }
    }
    let priorities: Vec<f64> = (0..n).map(|i| if i == ref_i { 1.0 } else { 0.0 }).collect();
    if best - second < JAVA_INFORMATIVE {
        let mut best_pri = priorities[best_i];
        let mut second_pri = priorities[second_i];
        for a in 0..n {
            let cand = lls[a];
            if a == best_i || best - cand > JAVA_INFORMATIVE {
                continue;
            }
            let pri = priorities[a];
            if pri > best_pri {
                second_i = best_i;
                best_i = a;
                second_pri = best_pri;
                best_pri = pri;
            } else if pri > second_pri {
                second_i = a;
                second_pri = pri;
            }
        }
    }
    let best_ll = lls[best_i];
    let second_ll = if second_i != best_i {
        lls[second_i]
    } else {
        f64::NEG_INFINITY
    };
    let conf = if best_ll == second_ll {
        0.0
    } else {
        best_ll - second_ll
    };
    (best_i, second_i, conf, conf > JAVA_INFORMATIVE)
}

fn simple_vote(lls: &[f64]) -> Option<usize> {
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
    format!("{}|{}|{}|{:?}", bases, bq, rec.flags(), rec.cigar())
}

fn fnv1a64(s: &str) -> u64 {
    let mut h = 0xcbf29ce484222325u64;
    for b in s.as_bytes() {
        h ^= u64::from(*b);
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

/// Identity remarg is column select (max of a singleton), not log-sum.
#[test]
fn forensic_6r89_remarg_is_max_of_mapped_columns_not_log_sum() {
    let old = [0.0, -3.0, -1.0, -8.0];
    let keep = remarg_identity(&old, &[0, 2]);
    assert_eq!(keep, vec![0.0, -1.0]);
    let log_sum = old[0] + (10.0_f64.powf(old[2] - old[0])).log10();
    assert_ne!(keep[0], log_sum);
}

/// DepthPerAlleleBySample maps each remaining allele to itself only (no pooling of unused ALTs).
#[test]
fn forensic_6r89_remaining_allele_map_is_identity_not_pool_unused() {
    // 0=TG, 1=T, 2=CG, 3=*
    let old = [-2.0, 0.0, -4.0, -1.0];
    let identity = remarg_identity(&old, &[0, 2]);
    assert_eq!(identity, vec![-2.0, -4.0]);
    let pooled_ref = old[0].max(old[1]).max(old[3]);
    assert_eq!(pooled_ref, 0.0);
    assert_ne!(identity[0], pooled_ref);
}

/// `isInformative` is strict `confidence > 0.2`. Equality is uninformative.
#[test]
fn forensic_6r89_is_informative_is_strict_greater_than_0_2() {
    let eq = [0.0, -0.2];
    let (_, _, conf, inf) = java_best_breaking_ties(&eq, 0);
    assert!((conf - 0.2).abs() < 1e-15);
    assert!(!inf, "confidence == 0.2 is not informative");
    let over = [0.0, -0.2000001];
    let (_, _, conf2, inf2) = java_best_breaking_ties(&over, 0);
    assert!(conf2 > 0.2);
    assert!(inf2);
}

/// Tie-break: when (best-second) < 0.2, REF priority wins assignment; confidence uses post-swap likelihoods.
#[test]
fn forensic_6r89_tie_break_prefers_ref_when_gap_below_threshold() {
    // ALT slightly better (gap 0.1 < 0.2) → assign REF; confidence negative → uninformative.
    let lls = [-0.1, 0.0];
    let (best, _, conf, inf) = java_best_breaking_ties(&lls, 0);
    assert_eq!(best, 0, "REF priority wins the near-tie");
    assert!(conf < 0.0);
    assert!(!inf);
    let simple = simple_vote(&lls);
    assert_eq!(simple, None, "simple vote also uninformative at gap 0.1");
}

/// Gap ≥ 0.2 does not enter the priority path; ALT keeps the win if gap > 0.2.
#[test]
fn forensic_6r89_tie_break_does_not_steal_clear_alt() {
    let lls = [-0.5, 0.0];
    let (best, _, conf, inf) = java_best_breaking_ties(&lls, 0);
    assert_eq!(best, 1);
    assert!((conf - 0.5).abs() < 1e-12);
    assert!(inf);
}

/// Extra 4-way allele floor before remarg can collapse remaining-allele gaps.
/// Java DepthPerAlleleBySample remargs the allele matrix as stored (haplotype-normalized only).
#[test]
fn forensic_6r89_four_way_floor_before_remarg_can_collapse_remaining_gap() {
    let mut four = vec![-9.0, 0.0, -10.0, -8.0]; // T wins 4-way; raw TG vs CG still 1.0 apart (REF)
    let raw = remarg_identity(&four, &[0, 2]);
    assert_eq!(simple_vote(&raw), Some(0));
    apply_allele_floor(&mut four);
    let floored = remarg_identity(&four, &[0, 2]);
    assert_eq!(four[0], -4.5);
    assert_eq!(four[2], -4.5);
    assert_eq!(simple_vote(&floored), None);
    assert_ne!(simple_vote(&raw), simple_vote(&floored));
}

/// Identity remarg of two remaining columns cannot assign more AD than the number of
/// rows whose remaining likelihoods are unequal. Java AD 36+19=55 is not this operation
/// on a 62-row TG-vs-CG matrix that contains 20 exact remaining-allele ties.
#[test]
fn forensic_6r89_identity_remarg_cannot_exceed_unequal_remaining_rows() {
    let rows = [
        [0.0, -1.0],
        [-1.0, 0.0],
        [0.0, 0.0],
        [0.0, 0.0],
        [0.0, -0.01],
    ];
    let mut unequal = 0i32;
    let mut ad = [0i32; 2];
    for two in &rows {
        if two[0] != two[1] {
            unequal += 1;
        }
        if let Some(v) = simple_vote(two) {
            ad[v] += 1;
        }
    }
    assert_eq!(unequal, 3);
    assert!(ad[0] + ad[1] <= unequal);
    assert_ne!(
        ad,
        [3, 2],
        "must not invent counts from exact remaining-allele ties"
    );
}
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

fn ad_from_two(rows: &[[f64; 2]], use_java_ties: bool) -> [i32; 2] {
    let mut ad = [0i32; 2];
    for two in rows {
        if use_java_ties {
            let (best, _, _, inf) = java_best_breaking_ties(two, 0);
            if inf {
                ad[best] += 1;
            }
        } else if let Some(v) = simple_vote(two) {
            ad[v] += 1;
        }
    }
    ad
}

#[test]
fn live_remaining_allele_remarg_and_informative() {
    if std::env::var("HOLDOUT_6R89").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R89=1");
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
    let snap = take_colocated_merge_numerics()
        .into_iter()
        .find(|n| n.loc == POS_SNP)
        .expect("snap");

    let haps = &outcome.assembly.haplotypes;
    let full_pad = outcome.assembly.padded_reference_start_1based();
    let full_ref = outcome.assembly.reference_bases();
    let hap_events = build_per_haplotype_variation_events(
        haps,
        full_ref,
        full_pad,
        outcome.assembly.max_mnp_distance(),
        covering[0].contig.as_str(),
    );
    let _replaced = replace_span_del_events(
        &variation_events_at_position_from_cache(&hap_events, POS_SNP, true),
        POS_SNP,
        full_pad,
        full_ref,
    );

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
    let cg_ai = alts.iter().position(|a| a == "CG").expect("CG");
    let t_ai = alts.iter().position(|a| a == "T");
    let star_ai = alts.iter().position(|a| a == SPAN_DEL_ALLELE);

    let loc = POS_SNP;
    let end = loc.saturating_add(long_ref.len().saturating_sub(1) as u64);
    let margin = DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN;
    let pairhmm_idx: std::collections::HashSet<usize> = outcome
        .read_likelihoods
        .iter()
        .map(|e| e.read_index.get())
        .collect();
    let overlap_ll: Vec<_> = outcome
        .read_likelihoods
        .iter()
        .filter(|e| {
            let i = e.read_index.get();
            pairhmm_idx.contains(&i)
                && outcome
                    .genotyping_reads
                    .get(i)
                    .is_some_and(|r| java_alignment_read_overlaps_interval(r, loc, end, margin))
        })
        .cloned()
        .collect();
    let hap_rows = region_likelihoods_to_rows(&overlap_ll, haps.len());
    let mut four_raw: Vec<(usize, Vec<f64>)> = hap_rows
        .iter()
        .map(|row| {
            let lls: Vec<f64> = pools.iter().map(|p| pool_max(row, p)).collect();
            (row.read_index, lls)
        })
        .collect();
    four_raw.sort_by_key(|(i, _)| *i);
    assert_eq!(four_raw.len(), 62);

    let mut two_raw = Vec::new();
    let mut two_floor = Vec::new();
    let mut two_pool_unused_into_ref = Vec::new();
    let mut n_floor_clips = 0usize;
    let mut n_tie_gap_lt_0_2 = 0usize;
    let mut n_java_tie_flips = 0usize;
    let mut informative_deltas: Vec<(f64, usize, bool)> = Vec::new();
    let mut uninf_deltas: Vec<(f64, usize)> = Vec::new();
    let mut diverge: Vec<String> = Vec::new();

    for (ri, raw) in &four_raw {
        let mut floored = raw.clone();
        let before = floored.clone();
        apply_allele_floor(&mut floored);
        if floored
            .iter()
            .zip(before.iter())
            .any(|(a, b)| (*a - *b).abs() > 0.0)
        {
            n_floor_clips += 1;
        }
        let raw_two = [raw[0], raw[cg_ai + 1]];
        let floor_two = [floored[0], floored[cg_ai + 1]];
        two_raw.push(raw_two);
        two_floor.push(floor_two);

        let mut pooled_ref = raw[0];
        if let Some(ti) = t_ai {
            pooled_ref = pooled_ref.max(raw[ti + 1]);
        }
        if let Some(si) = star_ai {
            pooled_ref = pooled_ref.max(raw[si + 1]);
        }
        two_pool_unused_into_ref.push([pooled_ref, raw[cg_ai + 1]]);

        let (jb, _, jconf, jinf) = java_best_breaking_ties(&raw_two, 0);
        let sv = simple_vote(&raw_two);
        if (raw_two[0] - raw_two[1]).abs() < JAVA_INFORMATIVE {
            n_tie_gap_lt_0_2 += 1;
        }
        if sv != jinf.then_some(jb) {
            n_java_tie_flips += 1;
        }
        if jinf {
            informative_deltas.push((jconf, *ri, jb == 0));
        } else {
            uninf_deltas.push((jconf.abs(), *ri));
        }
        let rust_floor_vote = simple_vote(&floor_two);
        if sv != rust_floor_vote || jinf.then_some(jb) != rust_floor_vote {
            let rec = &outcome.genotyping_reads[*ri];
            diverge.push(format!(
                "fp={:016x} rawTG={:.12} rawCG={:.12} floorTG={:.12} floorCG={:.12} java_best={} java_inf={} rust_floor={:?} simple_raw={:?}",
                fnv1a64(&read_fp_key(rec)),
                raw_two[0],
                raw_two[1],
                floor_two[0],
                floor_two[1],
                jb,
                jinf,
                rust_floor_vote,
                sv,
            ));
        }
    }

    let ad_raw_simple = ad_from_two(&two_raw, false);
    let ad_raw_java_ties = ad_from_two(&two_raw, true);
    let ad_floor_simple = ad_from_two(&two_floor, false);
    let ad_floor_java_ties = ad_from_two(&two_floor, true);
    let ad_pool_unused = ad_from_two(&two_pool_unused_into_ref, true);

    let mut ad_4way_any_map_unused_to_ref = [0i32; 2];
    let mut ad_4way_inf_map_unused_to_ref = [0i32; 2];
    for (_ri, raw) in &four_raw {
        let mut best_i = 0usize;
        let mut best = f64::NEG_INFINITY;
        let mut second = f64::NEG_INFINITY;
        for (i, &ll) in raw.iter().enumerate() {
            if ll > best {
                second = best;
                best = ll;
                best_i = i;
            } else if ll > second {
                second = ll;
            }
        }
        let mapped = if best_i == cg_ai + 1 { 1 } else { 0 };
        ad_4way_any_map_unused_to_ref[mapped] += 1;
        if best.is_finite() && (best - second).abs() > LOG_10_INFORMATIVE_THRESHOLD {
            ad_4way_inf_map_unused_to_ref[mapped] += 1;
        }
    }
    let mut ad_any_gap = [0i32; 2];
    for two in &two_raw {
        if two[0] > two[1] {
            ad_any_gap[0] += 1;
        } else if two[1] > two[0] {
            ad_any_gap[1] += 1;
        }
    }

    informative_deltas.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    uninf_deltas.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap());

    let vcf_ad = vcf.samples[0].ad.clone().unwrap_or_default();
    let vcf_pl = vcf.samples[0].pl.clone().unwrap_or_default();
    eprintln!(
        "6R.89 alleles_4way=[{:?}, {:?}] remaining=[TG, CG] n={} floor_clips={} snap_floor_clips={} ties_gap<0.2={} java_tie_assignment_flips={}",
        long_ref,
        alts,
        four_raw.len(),
        n_floor_clips,
        snap.n_allele_floor_clips,
        n_tie_gap_lt_0_2,
        n_java_tie_flips
    );
    eprintln!(
        "6R.89 AD raw_simple={:?} raw_java_ties={:?} floor_simple={:?} floor_java_ties={:?} pool_unused_into_ref={:?} any_strict_best={:?} four_any_map_unused_to_ref={:?} four_inf_map_unused_to_ref={:?} snap_remarg={:?} snap_perm={:?} vcf={:?} java_oracle=[36, 19]",
        ad_raw_simple,
        ad_raw_java_ties,
        ad_floor_simple,
        ad_floor_java_ties,
        ad_pool_unused,
        ad_any_gap,
        ad_4way_any_map_unused_to_ref,
        ad_4way_inf_map_unused_to_ref,
        snap.subset_ad_remarginalized,
        snap.subset_ad_permuted,
        vcf_ad
    );
    eprintln!("6R.89 diverge_raw_vs_floor n={}", diverge.len());
    for (i, line) in diverge.iter().take(16).enumerate() {
        eprintln!("6R.89 DIVERGE[{i}] {line}");
    }
    eprintln!("6R.89 smallest_informative_deltas:");
    for (d, ri, is_ref) in informative_deltas.iter().take(10) {
        eprintln!("  delta={d:.12} ri={ri} ref={is_ref}");
    }
    eprintln!("6R.89 largest_uninformative_deltas:");
    for (d, ri) in uninf_deltas.iter().take(10) {
        eprintln!("  |delta|={d:.12} ri={ri}");
    }
    eprintln!(
        "6R.89 vcf GT={:?} AD={:?} PL={:?} QUAL={:?} sample=NA12878 n_samples=1",
        vcf.samples[0].gt.as_ref().map(|g| g.alleles.clone()),
        vcf_ad,
        vcf_pl,
        vcf.quality
    );

    assert_eq!(snap.n_reads, 62);
    assert_eq!(ad_floor_simple.to_vec(), snap.subset_ad_remarginalized);
    assert_eq!(vcf_ad, vec![26u32, 9]);
    assert_eq!(
        ad_raw_simple, ad_raw_java_ties,
        "Test B: REF tie-break does not change informative AD on raw 2-way (near-ties stay uninformative)"
    );
    assert!(
        ad_any_gap[0] + ad_any_gap[1] < 55,
        "Java AD 36+19=55 exceeds the number of TG≠CG rows on equivalent C"
    );
}
