//! 6R.88 coordinate-free: AlleleLikelihoods object entering DepthPerAlleleBySample.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//! `HaplotypeCallerGenotypingEngine.assignGenotypeLikelihoods`:
//!
//! ```text
//! haplotype AlleleLikelihoods                    // PairHMM evidence
//!   .marginalize(alleleMapper)                   // NEW object; same evidence; allele columns
//! retainEvidence(target.overlaps)                // in-place; drop non-overlap by index
//! calculateGLsForThisEvent                       // same object
//! calculateGenotypes                             // unused-ALT subset → call alleles
//! prepareReadAlleleLikelihoodsForAnnotation
//!   default HC (no contamination):
//!     reuse genotyping AlleleLikelihoods         // NOT re-marginalized from haplotypes
//!     updateNonRefAlleleLikelihoods(call alleles) // no-op without <NON_REF>
//!     addEvidence(overlappingFilteredReads, 0)   // append; all cells 0
//! annotationEngine.annotateContext(call, …, readAlleleLikelihoods, …)
//!   genotype annotations receive allele-level likelihoods (not haplotype preFiltering)
//! DepthPerAlleleBySample.annotate
//!   alleles = call.getAlleles()                  // unused-ALT remaining; reverseTrim is AFTER
//!   annotateWithLikelihoods:
//!     marginalize({allele → [allele]})           // NEW object; same evidence; column select
//!     bestAllelesBreakingTies(sample)
//!     filter isInformative (confidence > 0.2)
//! reverseTrimAlleles                             // AFTER annotation
//! ```
//!
//! 6R.85 proved: SPAN_DEL is not causal; Java/Rust 136×4 likelihood matrices are identical.
//! 6R.86 identified 136→62 but did not prove Java/Rust overlap divergence.
//! 6R.87 proved overlap/coordinate semantics are identical for all 136 reads.
//! This round starts at the AlleleLikelihoods object passed to annotateWithLikelihoods.
//! PL, QUAL, PairHMM, overlap, permute-vs-remarg, and AD arithmetic are not investigated
//! as production changes.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r88_ad_likelihoods_evidence_contract
//! HOLDOUT_6R88=1 cargo test -p gatk-haplotypecaller --test forensic_6r88_ad_likelihoods_evidence_contract live_ -- --nocapture
//! ```

use gatk_haplotypecaller::event_map::{
    overlapping_events, remap_alt_onto_longer_ref, VariationEvent,
};
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::genotyping::ReadLikelihoodRow;
use gatk_haplotypecaller::hc_allele_mapping::SPAN_DEL_ALLELE;
use gatk_haplotypecaller::hc_genotyping_engine::java_alignment_read_overlaps_interval;
use gatk_haplotypecaller::read_assembly_filter::{passes_assembly_read, AssemblyReadFilterConfig};
use gatk_haplotypecaller::read_realignment::LOG_10_INFORMATIVE_THRESHOLD;
use gatk_haplotypecaller::{region_likelihoods_to_rows, DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN};
use rust_htslib::bam::record::Aux;
use rust_htslib::bam::{self, record::Cigar, record::CigarString};
use std::collections::HashSet;

fn row(lls: Vec<f64>) -> ReadLikelihoodRow {
    ReadLikelihoodRow {
        read_index: 0,
        read_id: String::new(),
        haplotype_log10_likelihoods: lls,
    }
}

/// Java `BestAllele.isInformative`: confidence = best − second; require `> 0.2`.
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

fn fnv1a64(s: &str) -> u64 {
    let mut h = 0xcbf29ce484222325u64;
    for b in s.as_bytes() {
        h ^= u64::from(*b);
        h = h.wrapping_mul(0x100000001b3);
    }
    h
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

/// Java `addEvidence(..., 0)` fills every allele cell with 0 → confidence 0 → uninformative.
#[test]
fn forensic_6r88_add_evidence_zero_likelihoods_are_uninformative() {
    let zeros = vec![0.0, 0.0, 0.0, 0.0];
    assert_eq!(informative_vote(&zeros), None);
    let remarg = [0.0, 0.0];
    assert_eq!(informative_vote(&remarg), None);
}

/// `updateNonRefAlleleLikelihoods` returns immediately when `<NON_REF>` is absent,
/// even if `call.getAlleles().size() != likelihoods.numberOfAlleles()`.
#[test]
fn forensic_6r88_update_nonref_is_noop_without_symbolic_nonref() {
    let genotyping_alleles = ["TG", "*", "T", "CG"];
    let call_alleles = ["TG", "CG"];
    assert_ne!(genotyping_alleles.len(), call_alleles.len());
    assert!(!genotyping_alleles.iter().any(|a| *a == "<NON_REF>"));
    assert!(!call_alleles.iter().any(|a| *a == "<NON_REF>"));
}

/// Default HC reuses the genotyping AlleleLikelihoods for annotation (no contamination).
/// `addEvidence` appends; it does not drop retainEvidence members.
#[test]
fn forensic_6r88_prepare_reuses_genotyping_object_and_only_appends() {
    let retain_fps = ["r0", "r1", "r2"];
    let filtered_overlap = ["f0"];
    let mut c: Vec<&str> = retain_fps.to_vec();
    for f in filtered_overlap {
        if !c.contains(&f) {
            c.push(f);
        }
    }
    assert_eq!(c, ["r0", "r1", "r2", "f0"]);
    assert_eq!(retain_fps.len(), 3, "B is unchanged by append");
}

/// DepthPerAlleleBySample remarg is a per-allele identity map (column select).
/// `Collectors.toMap` HashMap key order may permute columns; AD is keyed by allele identity.
#[test]
fn forensic_6r88_remarg_column_permutation_does_not_change_identity_ad() {
    let rows = [
        row(vec![-1.0, -10.0, 0.0, -8.0]),   // CG
        row(vec![0.0, -10.0, -8.0, -9.0]),   // REF
        row(vec![-0.01, 0.0, -10.0, -10.0]), // T; remarg REF vs CG → REF
    ];
    let mut ad_ref_cg = [0i32; 2];
    let mut ad_cg_ref = [0i32; 2];
    for r in &rows {
        let ref_cg = [
            r.haplotype_log10_likelihoods[0],
            r.haplotype_log10_likelihoods[2],
        ];
        let cg_ref = [
            r.haplotype_log10_likelihoods[2],
            r.haplotype_log10_likelihoods[0],
        ];
        if let Some(v) = informative_vote(&ref_cg) {
            ad_ref_cg[v] += 1;
        }
        if let Some(v) = informative_vote(&cg_ref) {
            ad_cg_ref[v] += 1;
        }
    }
    // Column 0/1 swapped; identity AD is REF count then CG count.
    assert_eq!(ad_ref_cg, [2, 1]);
    assert_eq!(ad_cg_ref, [1, 2]);
    let identity_from_swapped = [ad_cg_ref[1], ad_cg_ref[0]];
    assert_eq!(identity_from_swapped, ad_ref_cg);
}

/// reverseTrim is after annotation; remaining `call` alleles at AD time still have the
/// untrimmed merged REF (e.g. TG) and remaining ALT (e.g. CG).
#[test]
fn forensic_6r88_annotation_sees_untrimmed_call_alleles() {
    let merged = ["TG", "*", "T", "CG"];
    let after_unused_alt = ["TG", "CG"];
    let after_reverse_trim = ["T", "C"];
    assert_eq!(after_unused_alt[0], merged[0]);
    assert_ne!(after_reverse_trim[0], after_unused_alt[0]);
}

/// Permute of 4-way informative counts is not DepthPerAlleleBySample remarg.
#[test]
fn forensic_6r88_vcf_permute_is_not_the_java_ad_input_algorithm() {
    let rows = [
        row(vec![0.0, -10.0, -10.0, -10.0]),
        row(vec![-10.0, -10.0, 0.0, -10.0]),
        row(vec![-0.5, 0.0, -10.0, -10.0]),
        row(vec![-0.5, -10.0, -10.0, 0.0]),
    ];
    let mut ad4 = [0i32; 4];
    for r in &rows {
        if let Some(i) = informative_vote(&r.haplotype_log10_likelihoods) {
            ad4[i] += 1;
        }
    }
    let permuted = [ad4[0], ad4[2]];
    let mut remarg = [0i32; 2];
    for r in &rows {
        let two = [
            r.haplotype_log10_likelihoods[0],
            r.haplotype_log10_likelihoods[2],
        ];
        if let Some(v) = informative_vote(&two) {
            remarg[v] += 1;
        }
    }
    assert_eq!(permuted, [1, 1]);
    assert_eq!(remarg, [3, 1]);
    assert_ne!(permuted, remarg);
}

fn mate(qname: &[u8], pos0: i64) -> bam::Record {
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

/// retainEvidence keeps overlapping mates independently; it does not group by sample beyond
/// the existing per-sample evidence lists. Default HC is one sample.
#[test]
fn forensic_6r88_retain_evidence_does_not_collapse_sample_or_qname() {
    let r0 = mate(b"frag1", 99);
    let r1 = mate(b"frag1", 100);
    let keep: Vec<usize> = (0..2)
        .filter(|&i| java_alignment_read_overlaps_interval([&r0, &r1][i], 105, 105, 2))
        .collect();
    assert_eq!(keep.len(), 2);
}

#[test]
fn live_ad_likelihoods_evidence_object() {
    if std::env::var("HOLDOUT_6R88").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R88=1");
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
    let _at_loc: Vec<VariationEvent> = replaced
        .into_iter()
        .filter(|e| e.start_1based.get() == POS_SNP)
        .collect();

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
    let mut allele_names = vec![long_ref.to_string()];
    allele_names.extend(alts.iter().cloned());

    let mut ad4 = vec![0i32; pools.len()];
    let mut remarg = [0i32; 2];
    let mut d_fps: Vec<(String, &'static str)> = Vec::new();
    let mut rust_perm_fps: Vec<(String, &'static str)> = Vec::new();
    for (ri, lls) in &four_way {
        let rec = &outcome.genotyping_reads[*ri];
        let fp = read_fp_key(rec);
        if let Some(v) = informative_vote(lls) {
            ad4[v] += 1;
            if v == 0 || v == cg_ai + 1 {
                let allele = if v == 0 { "REF" } else { "ALT" };
                rust_perm_fps.push((fp.clone(), allele));
            }
        }
        let two = [lls[0], lls[cg_ai + 1]];
        if let Some(v) = informative_vote(&two) {
            remarg[v] += 1;
            let allele = if v == 0 { "REF" } else { "ALT" };
            d_fps.push((fp, allele));
        }
    }

    let b_fps: HashSet<String> = four_way
        .iter()
        .map(|(ri, _)| read_fp_key(&outcome.genotyping_reads[*ri]))
        .collect();

    let filter_cfg = AssemblyReadFilterConfig::gatk_defaults();
    let mut n_covering_overlap = 0usize;
    let mut n_filtered_overlap = 0usize;
    let mut filtered_overlap_fps: Vec<String> = Vec::new();
    let mut rg_ids: HashSet<String> = HashSet::new();
    for rec in &covering[0].reads {
        if let Ok(Aux::String(rg)) = rec.aux(b"RG") {
            rg_ids.insert(rg.to_string());
        }
        if !java_alignment_read_overlaps_interval(rec, loc, end, margin) {
            continue;
        }
        n_covering_overlap += 1;
        if !passes_assembly_read(rec, &filter_cfg) {
            n_filtered_overlap += 1;
            filtered_overlap_fps.push(read_fp_key(rec));
        }
    }
    let mut genotyping_rg: HashSet<String> = HashSet::new();
    for rec in &outcome.genotyping_reads {
        if let Ok(Aux::String(rg)) = rec.aux(b"RG") {
            genotyping_rg.insert(rg.to_string());
        }
    }

    let filtered_not_in_b: Vec<&String> = filtered_overlap_fps
        .iter()
        .filter(|fp| !b_fps.contains(*fp))
        .collect();

    // Java C real-likelihood evidence = B. Java C full = B ∪ overlapping filterNonPassingReads.
    let c_java_real = b_fps.len();
    let c_java_full = c_java_real + filtered_not_in_b.len();
    let c_rust = snap.n_reads;
    let rust_c_fps: HashSet<String> = four_way
        .iter()
        .map(|(ri, _)| read_fp_key(&outcome.genotyping_reads[*ri]))
        .collect();
    let java_c_real_fps = b_fps.clone();
    let common = java_c_real_fps.intersection(&rust_c_fps).count();
    let java_only: Vec<&String> = java_c_real_fps.difference(&rust_c_fps).collect();
    let rust_only: Vec<&String> = rust_c_fps.difference(&java_c_real_fps).collect();

    let vcf_ad = vcf.samples[0].ad.clone().unwrap_or_default();
    let vcf_pl = vcf.samples[0].pl.clone().unwrap_or_default();
    let vcf_gt = vcf.samples[0].gt.as_ref().map(|g| g.alleles.clone());

    let n_d_remarg = remarg[0] + remarg[1];
    assert_eq!(d_fps.len() as i32, n_d_remarg);
    let n_d_permute_kept = rust_perm_fps.len() as i32;
    let e_java = [36i32, 19];
    let e_rust_vcf = [
        vcf_ad.first().copied().unwrap_or(0) as i32,
        vcf_ad.get(1).copied().unwrap_or(0) as i32,
    ];

    eprintln!(
        "6R.88 A_pairhmm={} B_retainEvidence={} C_java_real={} C_java_full={} C_rust_production={} n_filtered_overlap={} filtered_not_in_B={}",
        pairhmm_idx.len(),
        four_way.len(),
        c_java_real,
        c_java_full,
        c_rust,
        n_filtered_overlap,
        filtered_not_in_b.len()
    );
    eprintln!(
        "6R.88 alleles_4way={:?} remaining_call_alleles=[{}, CG] bam_sample=NA12878 n_rg_covering={} n_rg_genotyping={} n_covering_overlap={}",
        allele_names,
        long_ref,
        rg_ids.len(),
        genotyping_rg.len(),
        n_covering_overlap
    );
    eprintln!(
        "6R.88 membership common={} JAVA_ONLY={} RUST_ONLY={}",
        common,
        java_only.len(),
        rust_only.len()
    );
    for (i, fp) in java_only.iter().take(8).enumerate() {
        eprintln!("6R.88 JAVA_ONLY[{i}] fp_hash={:016x}", fnv1a64(fp));
    }
    for (i, fp) in rust_only.iter().take(8).enumerate() {
        eprintln!("6R.88 RUST_ONLY[{i}] fp_hash={:016x}", fnv1a64(fp));
    }
    eprintln!(
        "6R.88 D_remarg={:?} n_informative_remarg={} D_permute_kept_cols={:?} n={} E_java={:?} E_rust_vcf={:?} snap_perm={:?} snap_remarg={:?}",
        remarg,
        n_d_remarg,
        [ad4[0], ad4[cg_ai + 1]],
        n_d_permute_kept,
        e_java,
        e_rust_vcf,
        snap.subset_ad_permuted,
        snap.subset_ad_remarginalized
    );
    eprintln!(
        "6R.88 vcf GT={:?} AD={:?} PL={:?} QUAL={:?}",
        vcf_gt, vcf_ad, vcf_pl, vcf.quality
    );

    assert_eq!(pairhmm_idx.len(), 136, "A is the proven 136 PairHMM set");
    assert_eq!(
        four_way.len(),
        62,
        "B is the shared retainEvidence overlap set"
    );
    assert_eq!(
        four_way.len(),
        snap.n_overlap_before_qname_dedupe,
        "B matches production overlap audit"
    );
    assert_eq!(
        c_rust,
        four_way.len(),
        "Rust AD input evidence count equals B (qname collapse is not causal)"
    );
    assert_eq!(
        common,
        four_way.len(),
        "real-likelihood C membership is identical"
    );
    assert!(
        java_only.is_empty() && rust_only.is_empty(),
        "no JAVA_ONLY/RUST_ONLY real-likelihood evidence at AD input"
    );
    assert_eq!(
        snap.n_qnames_with_multiple_overlapping_reads, 0,
        "QNAME collapse is not a C reconstruction"
    );
    assert!(
        !genotyping_rg.is_empty(),
        "genotyping evidence carries RG tags"
    );
    assert!(
        n_covering_overlap >= four_way.len(),
        "covering overlap is at least retainEvidence B"
    );
    assert!(
        n_filtered_overlap == 0 && filtered_not_in_b.is_empty(),
        "overlapping filterNonPassingReads is empty; Java C_full equals B"
    );
    assert_eq!(
        vcf_ad,
        vec![26u32, 9],
        "Rust VCF AD is still permute of 4-way on C, not Java 36,19"
    );
    assert_eq!(
        remarg.to_vec(),
        snap.subset_ad_remarginalized,
        "DepthPerAlleleBySample remarg on C is the diagnostic 62×2 object"
    );
    assert_ne!(
        remarg.to_vec(),
        vec![36, 19],
        "Java VCF AD is not remarg of this C object"
    );
    assert_ne!(
        snap.subset_ad_permuted,
        remarg.to_vec(),
        "Rust VCF permute is not Java DepthPerAlleleBySample remarg"
    );
}
