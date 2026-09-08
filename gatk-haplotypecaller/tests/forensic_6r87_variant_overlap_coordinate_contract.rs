//! 6R.87 coordinate-free: Java `variantCallingRelevantOverlap` / `target.overlaps(read)`.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//! `HaplotypeCallerGenotypingEngine.assignGenotypeLikelihoods`:
//!
//! ```text
//! mergedVC = makeMergedVariantContext(...)          // after replaceSpanDels
//! variantCallingRelevantOverlap =
//!     new SimpleInterval(mergedVC)
//!         .expandWithinContig(informativeReadOverlapMargin, dictionary)
//! // SimpleInterval is 1-based inclusive [start, end]
//! // expandWithinContig: start-padding, end+padding, clamp to contig
//! // default HC (not BQD/FRD):
//! //   (read, target) -> target.overlaps(read)
//! // SimpleInterval.overlaps (margin=0):
//! //   contig equal && this.start <= other.getEnd() && other.getStart() <= this.end
//! // GATKRead.getStart/getEnd = SAMRecord getAlignmentStart/getAlignmentEnd
//! // (1-based inclusive clipped alignment)
//! // Evidence objects are post-changeEvidence realigned reads.
//! ```
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r87_variant_overlap_coordinate_contract
//! HOLDOUT_6R87=1 cargo test -p gatk-haplotypecaller --test forensic_6r87_variant_overlap_coordinate_contract live_ -- --nocapture
//! ```

use gatk_haplotypecaller::hc_genotyping_engine::java_alignment_read_overlaps_interval;
use gatk_haplotypecaller::read_unclip::{alignment_end_1based, gatk_soft_start_1based};
use rust_htslib::bam::{self, record::Cigar, record::CigarString};
use std::collections::HashMap;

/// Java `SimpleInterval.expandWithinContig` then `overlaps` with margin 0.
/// Target and read coordinates are 1-based inclusive.
fn java_expanded_target_overlaps_read(
    merged_start: i64,
    merged_end: i64,
    padding: i32,
    contig_len: i64,
    read_start: i64,
    read_end: i64,
) -> bool {
    let (tstart, tend) = java_expand_within_contig(merged_start, merged_end, padding, contig_len);
    java_simple_interval_overlaps(tstart, tend, read_start, read_end)
}

fn java_expand_within_contig(start: i64, end: i64, padding: i32, contig_len: i64) -> (i64, i64) {
    let p = i64::from(padding.max(0));
    let bounded_start = (start - p).max(1);
    let bounded_stop = (end + p).min(contig_len);
    (bounded_start, bounded_stop)
}

/// Java `SimpleInterval.overlaps(Locatable)` with margin 0 (same contig assumed).
fn java_simple_interval_overlaps(tstart: i64, tend: i64, rstart: i64, rend: i64) -> bool {
    tstart <= rend && rstart <= tend
}

fn rec_align_1based(rec: &bam::Record) -> (i64, i64) {
    (
        (rec.pos() + 1).max(1),
        i64::from(alignment_end_1based(rec).max(1)),
    )
}

fn rec_soft_end_1based(rec: &bam::Record) -> i64 {
    let mut end = i64::from(alignment_end_1based(rec));
    for c in rec.cigar().iter().rev() {
        match c {
            Cigar::SoftClip(n) => end += i64::from(*n),
            Cigar::HardClip(_) => {}
            _ => break,
        }
    }
    end
}

fn mate(qname: &[u8], pos0: i64, cigar: CigarString) -> bam::Record {
    let mut rec = bam::Record::new();
    rec.set(qname, Some(&cigar), b"ACGTACGTAC", b"##########");
    rec.set_pos(pos0);
    rec
}

#[test]
fn forensic_6r87_simple_interval_is_1based_inclusive() {
    // [100, 102] size 3. Expand ±2 → [98, 104].
    let (s, e) = java_expand_within_contig(100, 102, 2, 1_000_000);
    assert_eq!((s, e), (98, 104));
}

#[test]
fn forensic_6r87_expand_is_symmetric_not_asymmetric() {
    let (s, e) = java_expand_within_contig(100, 101, 2, 1_000_000);
    assert_eq!(s, 100 - 2);
    assert_eq!(e, 101 + 2);
    let (s5, e5) = java_expand_within_contig(100, 101, 5, 1_000_000);
    assert_ne!((s, e), (s5, e5), "must not silently use ±5");
}

#[test]
fn forensic_6r87_overlaps_uses_alignment_not_softclip_for_default_hc() {
    // 10M ending 5bp before expanded start: alignment rejects; 5S trailing would accept if unclipped.
    let rec = mate(
        b"r",
        90,
        CigarString(vec![Cigar::Match(8), Cigar::SoftClip(5)]),
    );
    let (rs, re) = rec_align_1based(&rec);
    assert_eq!((rs, re), (91, 98));
    assert!(!java_simple_interval_overlaps(100, 104, rs, re));
    let soft_end = rec_soft_end_1based(&rec);
    assert_eq!(soft_end, 103);
    assert!(
        java_simple_interval_overlaps(100, 104, gatk_soft_start_1based(&rec).max(1), soft_end),
        "soft-unclip is BQD/FRD only; default HC must not use it"
    );
}

#[test]
fn forensic_6r87_zero_ref_consuming_read_end_equals_start() {
    let rec = mate(b"z", 99, CigarString(vec![Cigar::Ins(10)]));
    let (rs, re) = rec_align_1based(&rec);
    assert_eq!(rs, re);
    assert!(java_simple_interval_overlaps(100, 104, rs, re));
}

#[test]
fn forensic_6r87_rust_alignment_overlap_matches_java_expanded_target() {
    let rec = mate(b"r", 99, CigarString(vec![Cigar::Match(10)]));
    let merged_start = 105i64;
    let merged_end = 105i64;
    let padding = 2;
    let (rs, re) = rec_align_1based(&rec);
    let java =
        java_expanded_target_overlaps_read(merged_start, merged_end, padding, 1_000_000, rs, re);
    let rust = java_alignment_read_overlaps_interval(
        &rec,
        merged_start as u64,
        merged_end as u64,
        padding,
    );
    assert_eq!(java, rust);
    assert!(java, "10M at pos0=99 overlaps 105±2");
}

#[test]
fn forensic_6r87_gap5_before_expanded_start_is_reject_on_same_coords() {
    // Alignment end 5bp before expanded target start: same coords + Java predicate reject.
    // Do not widen padding to close a 5bp gap.
    let rec = mate(b"gap5", 83, CigarString(vec![Cigar::Match(10)]));
    let (rs, re) = rec_align_1based(&rec);
    assert_eq!((rs, re), (84, 93));
    let java = java_expanded_target_overlaps_read(100, 101, 2, 1_000_000, rs, re);
    let rust = java_alignment_read_overlaps_interval(&rec, 100, 101, 2);
    assert_eq!(java, rust);
    assert!(
        !java,
        "same coords + same predicate reject; do not widen padding"
    );
}

fn fp_key(rec: &bam::Record) -> String {
    let bases = String::from_utf8_lossy(&rec.seq().as_bytes()).into_owned();
    format!("{}|{}|{}", rec.qname().len(), rec.flags(), bases)
}

#[test]
fn live_variant_overlap_coordinates() {
    if std::env::var("HOLDOUT_6R87").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R87=1");
        return;
    }
    use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
    use gatk_haplotypecaller::event_map::{
        overlapping_events, remap_alt_onto_longer_ref, VariationEvent,
    };
    use gatk_haplotypecaller::genome_loc::GenomePosition;
    use gatk_haplotypecaller::hc_allele_mapping::SPAN_DEL_ALLELE;
    use gatk_haplotypecaller::read_realignment::LOG_10_INFORMATIVE_THRESHOLD;
    use gatk_haplotypecaller::{
        call_disposition, flatten_assembly_regions, region_likelihoods_to_rows,
        take_colocated_merge_numerics, traverse_assembly_region_walker,
        AssemblyRegionCallDisposition, CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams,
        WalkerTraversalConfig, DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN,
    };
    use std::collections::HashSet;
    use std::path::Path;

    const INTERVAL: &str = "20:29455000-29456500";
    const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
    const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
    const POS_SNP: u64 = 29_456_344;
    const CONTIG_LEN: i64 = 63_025_520;

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
    let orig_by_fp: HashMap<String, (i64, i64, i64, i64, String)> = covering[0]
        .reads
        .iter()
        .map(|r| {
            let (s, e) = rec_align_1based(r);
            (
                fp_key(r),
                (
                    s,
                    e,
                    gatk_soft_start_1based(r).max(1),
                    rec_soft_end_1based(r),
                    format!("{:?}", r.cigar()),
                ),
            )
        })
        .collect();

    let outcome = HaplotypeCallerEngine::call_region(
        covering[0],
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let snap = take_colocated_merge_numerics()
        .into_iter()
        .find(|n| n.loc == POS_SNP)
        .expect("colocated merge numerics");

    let merged_start = POS_SNP as i64;
    let merged_end = merged_start + snap.long_ref.len() as i64 - 1;
    let padding = DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN;
    let (tstart, tend) = java_expand_within_contig(merged_start, merged_end, padding, CONTIG_LEN);
    eprintln!(
        "6R.87 JAVA_TARGET contig=20 mergedVC={}-{} long_ref={} expand±{} → {}-{} (1-based inclusive)",
        merged_start, merged_end, snap.long_ref, padding, tstart, tend
    );

    let pairhmm_idx: HashSet<usize> = outcome
        .read_likelihoods
        .iter()
        .map(|e| e.read_index.get())
        .collect();
    let haps = &outcome.assembly.haplotypes;
    let loc = POS_SNP;
    let end = loc.saturating_add(snap.long_ref.len().saturating_sub(1) as u64);

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
        if spanning.is_empty() {
            return JavaPool::Ref;
        }
        let loc_pos = GenomePosition::new_1based(loc);
        for ev in spanning {
            if ev.start_1based == loc_pos {
                let alt = if ev.ref_allele == long_ref {
                    ev.alt_allele.clone()
                } else {
                    remap_alt_onto_longer_ref(&ev.ref_allele, &ev.alt_allele, long_ref)
                        .unwrap_or_else(|| ev.alt_allele.clone())
                };
                if let Some(ai) = alts.iter().position(|a| a == &alt) {
                    return JavaPool::Alt(ai);
                }
            } else if ev.start_1based < loc_pos {
                return if emit_spanning {
                    JavaPool::SpanDel
                } else {
                    JavaPool::Ref
                };
            }
        }
        JavaPool::Unassigned
    }
    fn pool_max(row: &gatk_haplotypecaller::genotyping::ReadLikelihoodRow, pool: &[usize]) -> f64 {
        pool.iter()
            .filter_map(|&i| row.haplotype_log10_likelihoods.get(i).copied())
            .fold(f64::NEG_INFINITY, f64::max)
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

    let alts = snap.alts.as_slice();
    let long_ref = snap.long_ref.as_str();
    let hap_events = gatk_haplotypecaller::event_map::build_per_haplotype_variation_events(
        haps,
        outcome.assembly.reference_bases(),
        outcome.assembly.padded_reference_start_1based(),
        outcome.assembly.max_mnp_distance(),
        covering[0].contig.as_str(),
    );
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
    let cg_ai = alts.iter().position(|a| a == "CG").unwrap_or(0);

    let mut n_rust_overlap = 0usize;
    let mut n_java_realigned = 0usize;
    let mut n_java_original = 0usize;
    let mut n_java_soft_realigned = 0usize;
    let mut first_java_true_rust_false: Option<usize> = None;
    let mut first_rust_true_java_false: Option<usize> = None;
    let mut orig_overlap_rust_reject_informative: Vec<usize> = Vec::new();
    let mut idx59: Option<String> = None;

    let all_ll: Vec<_> = outcome
        .read_likelihoods
        .iter()
        .filter(|e| pairhmm_idx.contains(&e.read_index.get()))
        .cloned()
        .collect();
    let all_rows = region_likelihoods_to_rows(&all_ll, haps.len());
    let remarg_of = |row: &gatk_haplotypecaller::genotyping::ReadLikelihoodRow| -> Option<usize> {
        let mut lls: Vec<f64> = pools.iter().map(|p| pool_max(row, p)).collect();
        if let Some(&best) = lls
            .iter()
            .filter(|v| v.is_finite())
            .max_by(|a, b| a.total_cmp(b))
        {
            let floor = best - 4.5;
            for v in &mut lls {
                if v.is_finite() && *v < floor {
                    *v = floor;
                }
            }
        }
        informative_vote(&[
            lls[0],
            lls.get(cg_ai + 1).copied().unwrap_or(f64::NEG_INFINITY),
        ])
    };

    for &i in &{
        let mut v: Vec<usize> = pairhmm_idx.iter().copied().collect();
        v.sort_unstable();
        v
    } {
        let Some(rec) = outcome.genotyping_reads.get(i) else {
            continue;
        };
        let (rs, re) = rec_align_1based(rec);
        let rust_ov = java_alignment_read_overlaps_interval(rec, loc, end, padding);
        let java_realigned = java_expanded_target_overlaps_read(
            merged_start,
            merged_end,
            padding,
            CONTIG_LEN,
            rs,
            re,
        );
        let java_soft = java_expanded_target_overlaps_read(
            merged_start,
            merged_end,
            padding,
            CONTIG_LEN,
            gatk_soft_start_1based(rec).max(1),
            rec_soft_end_1based(rec),
        );
        if rust_ov {
            n_rust_overlap += 1;
        }
        if java_realigned {
            n_java_realigned += 1;
        }
        if java_soft {
            n_java_soft_realigned += 1;
        }
        if java_realigned && !rust_ov && first_java_true_rust_false.is_none() {
            first_java_true_rust_false = Some(i);
        }
        if rust_ov && !java_realigned && first_rust_true_java_false.is_none() {
            first_rust_true_java_false = Some(i);
        }

        let orig = orig_by_fp.get(&fp_key(rec)).cloned();
        let (os, oe, osoft_s, osoft_e, ocigar) = orig.unwrap_or((rs, re, rs, re, String::new()));
        let java_original = java_expanded_target_overlaps_read(
            merged_start,
            merged_end,
            padding,
            CONTIG_LEN,
            os,
            oe,
        );
        if java_original {
            n_java_original += 1;
        }

        let row = all_rows.iter().find(|r| r.read_index == i);
        let remarg = row.and_then(remarg_of);
        if java_original && !rust_ov && remarg.is_some() {
            orig_overlap_rust_reject_informative.push(i);
        }

        if i == 59 {
            idx59 = Some(format!(
                "idx=59 rust_overlap={rust_ov} java_realigned={java_realigned} java_original={java_original} java_soft_realigned={java_soft} \
                 realigned={rs}-{re} original={os}-{oe} orig_soft={osoft_s}-{osoft_e} realigned_soft={}-{} \
                 realigned_cigar={:?} orig_cigar={ocigar} remarg={remarg:?} qname_len={} flags={}",
                gatk_soft_start_1based(rec).max(1),
                rec_soft_end_1based(rec),
                rec.cigar(),
                rec.qname().len(),
                rec.flags()
            ));
        }
    }

    eprintln!("6R.87 READ59 {}", idx59.as_deref().unwrap_or("missing"));
    eprintln!(
        "6R.87 COUNTS n_pairhmm={} rust_overlap={} java_on_realigned={} java_on_original={} java_soft_realigned={}",
        pairhmm_idx.len(),
        n_rust_overlap,
        n_java_realigned,
        n_java_original,
        n_java_soft_realigned
    );
    eprintln!(
        "6R.87 FIRST_JAVA_TRUE_RUST_FALSE={:?} FIRST_RUST_TRUE_JAVA_FALSE={:?} orig_overlap_but_rust_reject_informative={:?}",
        first_java_true_rust_false,
        first_rust_true_java_false,
        orig_overlap_rust_reject_informative
    );

    assert_eq!(
        n_rust_overlap, n_java_realigned,
        "Rust alignment overlap must match Java SimpleInterval(mergedVC).expand(2).overlaps(realigned read)"
    );
    assert!(
        first_java_true_rust_false.is_none() && first_rust_true_java_false.is_none(),
        "no Java/Rust overlap disagreement on post-realign genotyping reads"
    );
}
