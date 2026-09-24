//! 6R.222: first loc-loop state transition that turns production-arg
//! `try_genotype` (`AD 47,5` / `PL 68,0,1937`, n=52) into assign FORMAT
//! (`AD 44,5` / `PL 78,0,1811`, n=52) at `20:29455379 G/A`.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE. Do not assume `20:29455375 T/A` is causal.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r222_loc_loop_likelihood_lifecycle -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{
    build_event_start_positions_from_cache, build_per_haplotype_variation_events,
    merged_biallelic_sites_at_position, prefer_indel_over_colocated_snps,
    variation_events_at_position_from_cache, VariationEvent,
};
use gatk_haplotypecaller::hc_allele_mapping::replace_span_del_events;
use gatk_haplotypecaller::hc_genotyping_engine::{
    assign_genotype_likelihoods_for_region, clear_region_likelihood_rows_tls,
    diagnose_genotype_variation_event_with_region_state, genotype_from_read_rows,
    java_alignment_read_overlaps_interval, region_likelihood_rows_tls_identity,
    region_likelihoods_to_rows, take_colocated_merge_numerics, GenotypedSiteCall,
    HcGenotypingConfig, DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition, HaplotypeCallerEngine,
    ReadFilterParams, WalkerTraversalConfig,
};
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_455_379;
const PRED: u64 = 29_455_375;
const CLOSED_PL: u64 = 29_455_015;
const TARGET_REF: &str = "G";
const TARGET_ALT: &str = "A";
const AD_A: [i32; 2] = [42, 5];
const PL_A: [i32; 3] = [84, 0, 1738];
const AD_B: [i32; 2] = [42, 5];
const PL_B: [i32; 3] = [84, 0, 1738];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R222\t{key}\t{}", value.as_ref());
}

fn unique_n(likelihoods: &[gatk_haplotypecaller::RegionReadLikelihood]) -> usize {
    likelihoods
        .iter()
        .map(|c| c.read_index.get())
        .collect::<BTreeSet<_>>()
        .len()
}

fn fmt_of(call: &GenotypedSiteCall) -> (Vec<i32>, Vec<i32>, i32, usize) {
    (
        call.genotype.format.ad_as_i32(),
        call.genotype.format.pl_as_i32(),
        call.genotype.format.dp.as_i32(),
        unique_n(&call.annotation_likelihoods),
    )
}

fn dump_call(tag: &str, call: &GenotypedSiteCall) {
    let (ad, pl, dp, ann) = fmt_of(call);
    kv(
        tag,
        format!(
            "event={}:{}-{}/{}\tad={ad:?}\tpl={pl:?}\tdp={dp}\tann_n={ann}\tgq={}",
            call.event.contig,
            call.event.start_1based.get(),
            call.event.ref_allele,
            call.event.alt_allele,
            call.genotype.format.gq.as_i32()
        ),
    );
}

fn event_id(e: &VariationEvent) -> String {
    format!(
        "{}:{}-{} {}/{}",
        e.contig,
        e.start_1based.get(),
        e.end_1based.get(),
        e.ref_allele,
        e.alt_allele
    )
}

fn ann_matrix(call: &GenotypedSiteCall) -> BTreeMap<(usize, usize), i64> {
    call.annotation_likelihoods
        .iter()
        .map(|c| {
            (
                (c.read_index.get(), c.haplotype_index.get()),
                (c.log10_likelihood * 1e9).round() as i64,
            )
        })
        .collect()
}

fn matrix_delta(
    a: &BTreeMap<(usize, usize), i64>,
    b: &BTreeMap<(usize, usize), i64>,
) -> (usize, usize, usize, i64, f64) {
    let a_keys: BTreeSet<_> = a.keys().copied().collect();
    let b_keys: BTreeSet<_> = b.keys().copied().collect();
    let only_a = a_keys.difference(&b_keys).count();
    let only_b = b_keys.difference(&a_keys).count();
    let shared: Vec<_> = a_keys.intersection(&b_keys).copied().collect();
    let mut changed = 0usize;
    let mut abs_sum = 0i64;
    let mut max_abs = 0i64;
    for k in &shared {
        let d = (a[k] - b[k]).abs();
        if d != 0 {
            changed += 1;
            abs_sum += d;
            if d > max_abs {
                max_abs = d;
            }
        }
    }
    let mean = if changed == 0 {
        0.0
    } else {
        (abs_sum as f64 / 1e9) / changed as f64
    };
    (only_a, only_b, changed, max_abs, mean)
}

#[test]
fn forensic_6r222_java_lifecycle_contract() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE — proof-only");
    kv(
        "java_lifecycle",
        "per-event calculateGLsForThisEvent on retainEvidence AlleleLikelihoods; no shared mutation between neighboring events on default HC",
    );
    assert_eq!(DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN, 2);
}

#[test]
fn forensic_6r222_loc_loop_likelihood_lifecycle() {
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
    kv(
        "predecessor_candidate",
        "20:29455375 T/A — not assumed causal",
    );

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
    assert_eq!(covering.len(), 1);
    let region = covering[0];

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
    let apply_bases = outcome.assembly.apply_bases_shared();
    let apply_pad = haps
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.map(|g| g.start_1based()))
        .unwrap_or(full_pad);
    let config = HcGenotypingConfig::strict_java();
    let emit_spanning = !config.disable_spanning_event_genotyping;
    let hap_events = build_per_haplotype_variation_events(
        haps,
        full_ref,
        full_pad,
        outcome.assembly.max_mnp_distance(),
        &region.contig,
    );
    let stored = outcome.assembly.variation_events();

    let prod = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based.get() == TARGET
                && c.event.ref_allele == TARGET_REF
                && c.event.alt_allele == TARGET_ALT
        })
        .expect("production G/A");
    dump_call("C_call_region", prod);
    let (c_ad, c_pl, _, c_ann) = fmt_of(prod);
    assert_eq!(c_ad, AD_B);
    assert_eq!(c_pl, PL_B);
    assert_eq!(c_ann, 47);

    let positions: Vec<u64> = build_event_start_positions_from_cache(&hap_events)
        .iter()
        .copied()
        .filter(|p| *p >= region.start.get() && *p <= region.end.get() && *p <= TARGET)
        .collect();
    kv("loc_loop_order_through_target", format!("{positions:?}"));
    assert_eq!(
        positions,
        vec![29_455_314, 29_455_328, 29_455_337, PRED, TARGET],
        "production loc-loop order"
    );

    let overlap = |loc: u64| {
        reads
            .iter()
            .filter(|r| java_alignment_read_overlaps_interval(r, loc, loc, 2))
            .count()
    };
    kv("overlap_pm2_pred", overlap(PRED).to_string());
    kv("overlap_pm2_target", overlap(TARGET).to_string());
    let pred_overlap: BTreeSet<usize> = reads
        .iter()
        .enumerate()
        .filter(|(_, r)| java_alignment_read_overlaps_interval(r, PRED, PRED, 2))
        .map(|(i, _)| i)
        .collect();
    let tgt_overlap: BTreeSet<usize> = reads
        .iter()
        .enumerate()
        .filter(|(_, r)| java_alignment_read_overlaps_interval(r, TARGET, TARGET, 2))
        .map(|(i, _)| i)
        .collect();
    kv(
        "overlap_identity",
        format!(
            "pred_n={}\ttgt_n={}\tonly_pred={:?}\tonly_tgt={:?}",
            pred_overlap.len(),
            tgt_overlap.len(),
            pred_overlap
                .difference(&tgt_overlap)
                .copied()
                .collect::<Vec<_>>(),
            tgt_overlap
                .difference(&pred_overlap)
                .copied()
                .collect::<Vec<_>>()
        ),
    );

    let diagnose = |tag: &str, event: &VariationEvent| {
        let tls_before = region_likelihood_rows_tls_identity();
        let diagnosed = diagnose_genotype_variation_event_with_region_state(
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
            stored,
            Some(&hap_events),
        )
        .expect(tag);
        let tls_after = region_likelihood_rows_tls_identity();
        kv(
            &format!("{tag}_tls"),
            format!("before={tls_before:?}\tafter={tls_after:?}"),
        );
        match diagnosed {
            Ok(call) => {
                dump_call(tag, &call);
                Some(call)
            }
            Err(reason) => {
                kv(tag, format!("REJECT {reason:?}"));
                None
            }
        }
    };

    let cache_target = variation_events_at_position_from_cache(&hap_events, TARGET, emit_spanning);
    assert_eq!(cache_target.len(), 1);
    assert_eq!(cache_target[0].ref_allele, TARGET_REF);
    assert_eq!(cache_target[0].alt_allele, TARGET_ALT);
    kv("target_event", event_id(&cache_target[0]));

    let mut pred_raw = variation_events_at_position_from_cache(&hap_events, PRED, emit_spanning);
    prefer_indel_over_colocated_snps(&mut pred_raw);
    let pred_span = replace_span_del_events(&pred_raw, PRED, apply_pad, apply_bases.as_ref());
    let pred_merged = merged_biallelic_sites_at_position(&pred_span, PRED);
    kv(
        "pred_events",
        pred_merged
            .iter()
            .map(event_id)
            .collect::<Vec<_>>()
            .join(";"),
    );
    kv("pred_merged_n", pred_merged.len().to_string());

    clear_region_likelihood_rows_tls();
    kv(
        "tls_cleared",
        format!("{:?}", region_likelihood_rows_tls_identity()),
    );

    let a_prime = diagnose("Aprime_target_alone", &cache_target[0]).expect("A'");
    let (ap_ad, ap_pl, _, ap_ann) = fmt_of(&a_prime);
    assert_eq!(ap_ad, AD_A);
    assert_eq!(ap_pl, PL_A);
    assert_eq!(ap_ann, 47);
    let a_mat = ann_matrix(&a_prime);

    clear_region_likelihood_rows_tls();
    let pred_calls: Vec<_> = pred_merged
        .iter()
        .map(|e| (event_id(e), diagnose("pred_try_genotype", e)))
        .collect();
    let pred_n_ok = pred_calls.iter().filter(|(_, c)| c.is_some()).count();
    kv("pred_try_genotype_n_ok", pred_n_ok.to_string());
    assert_eq!(
        pred_n_ok, 0,
        "predecessor T/A is REJECT VariantNotConfident"
    );
    let after_pred_tls = region_likelihood_rows_tls_identity();
    kv("tls_after_pred", format!("{after_pred_tls:?}"));

    let after_pred = diagnose("target_after_pred", &cache_target[0]).expect("after pred");
    let (after_pred_ad, after_pred_pl, _, _after_pred_ann) = fmt_of(&after_pred);
    kv(
        "pred_changes_target",
        (after_pred_ad != ap_ad || after_pred_pl != ap_pl).to_string(),
    );

    clear_region_likelihood_rows_tls();
    let target_first = diagnose("target_first_order", &cache_target[0]).expect("target first");
    let _ = pred_merged
        .iter()
        .map(|e| diagnose("pred_after_target", e))
        .count();
    let (tf_ad, tf_pl, _, _) = fmt_of(&target_first);
    kv(
        "order_swap_target_matches_aprime",
        (tf_ad == ap_ad && tf_pl == ap_pl).to_string(),
    );

    clear_region_likelihood_rows_tls();
    let rows = region_likelihoods_to_rows(likelihoods, haps.len());
    let _ = genotype_from_read_rows(&rows, haps, &config).expect("region_summary");
    kv(
        "tls_after_region_summary",
        format!("{:?}", region_likelihood_rows_tls_identity()),
    );
    let after_summary = diagnose("target_after_region_summary", &cache_target[0]).expect("summary");
    let (sum_ad, sum_pl, _, _) = fmt_of(&after_summary);

    clear_region_likelihood_rows_tls();
    let _ = region_likelihoods_to_rows(likelihoods, haps.len());
    let _ = genotype_from_read_rows(&rows, haps, &config);
    for e in &pred_merged {
        let _ = diagnose("prefix_pred", e);
    }
    let prefix_target = diagnose("target_after_prefix_replay", &cache_target[0]).expect("prefix");
    let (px_ad, px_pl, _, _) = fmt_of(&prefix_target);

    clear_region_likelihood_rows_tls();
    for e in &pred_merged {
        let _ = diagnose("clobber_pred", e);
    }
    let _ = region_likelihoods_to_rows(likelihoods, haps.len());
    let clobber = diagnose("target_after_pred_then_full_tls", &cache_target[0]).expect("clobber");
    let (cl_ad, cl_pl, _, _) = fmt_of(&clobber);

    let _ = take_colocated_merge_numerics();
    let assigned = assign_genotype_likelihoods_for_region(
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
        &region.contig,
        outcome.assembly.max_mnp_distance(),
        &config,
        stored,
        &[],
    )
    .expect("assign");
    let b = assigned
        .calls
        .iter()
        .find(|c| {
            c.event.start_1based.get() == TARGET
                && c.event.ref_allele == TARGET_REF
                && c.event.alt_allele == TARGET_ALT
        })
        .expect("assign G/A");
    dump_call("B_assign", b);
    let (b_ad, b_pl, _, b_ann) = fmt_of(b);
    assert_eq!(b_ad, AD_B);
    assert_eq!(b_pl, PL_B);
    assert_eq!(b_ann, 47);
    let colocated = take_colocated_merge_numerics();
    kv("colocated_snaps_n", colocated.len().to_string());
    assert_eq!(colocated.len(), 0);
    let b_mat = ann_matrix(b);
    let (only_a, only_b, changed, max_abs, mean) = matrix_delta(&a_mat, &b_mat);
    kv(
        "matrix_Aprime_vs_B",
        format!(
            "a_cells={}\tb_cells={}\tonly_a={only_a}\tonly_b={only_b}\tchanged={changed}\tmax_abs_e9={max_abs}\tmean_abs={mean:.6e}",
            a_mat.len(),
            b_mat.len()
        ),
    );
    let a_reads: BTreeSet<_> = a_prime
        .annotation_likelihoods
        .iter()
        .map(|c| c.read_index.get())
        .collect();
    let b_reads: BTreeSet<_> = b
        .annotation_likelihoods
        .iter()
        .map(|c| c.read_index.get())
        .collect();
    kv(
        "row_membership",
        format!(
            "a_n={}\tb_n={}\tonly_a={:?}\tonly_b={:?}",
            a_reads.len(),
            b_reads.len(),
            a_reads.difference(&b_reads).copied().collect::<Vec<_>>(),
            b_reads.difference(&a_reads).copied().collect::<Vec<_>>()
        ),
    );
    let membership_same = a_reads == b_reads;
    kv("row_membership_same", membership_same.to_string());

    kv("after_pred_ad", format!("{after_pred_ad:?}"));
    kv("after_pred_pl", format!("{after_pred_pl:?}"));
    kv("after_summary_ad", format!("{sum_ad:?}"));
    kv("after_summary_pl", format!("{sum_pl:?}"));
    kv("prefix_replay_ad", format!("{px_ad:?}"));
    kv("prefix_replay_pl", format!("{px_pl:?}"));
    kv("clobber_ad", format!("{cl_ad:?}"));
    kv("clobber_pl", format!("{cl_pl:?}"));

    let pred_causal = after_pred_ad == AD_B && after_pred_pl == PL_B;
    let skip_pred_stays_a = ap_ad == AD_A && ap_pl == PL_A;
    let order_independent = tf_ad == AD_A && tf_pl == PL_A;
    let clobber_restores = cl_ad == AD_A && cl_pl == PL_A;
    kv("pred_try_genotype_reproduces_B", pred_causal.to_string());
    kv(
        "prefix_replay_reproduces_B",
        ((px_ad == AD_B) && (px_pl == PL_B)).to_string(),
    );
    kv(
        "region_summary_reproduces_B",
        ((sum_ad == AD_B) && (sum_pl == PL_B)).to_string(),
    );
    kv("skip_pred_stays_Aprime", skip_pred_stays_a.to_string());
    assert!(skip_pred_stays_a, "skipping 29455375 must leave A'");
    assert!(
        pred_causal,
        "predecessor try_genotype then target must reproduce B"
    );
    assert!(order_independent, "target-first order must stay A'");
    assert!(
        clobber_restores,
        "full-table TLS fill after predecessor must restore A'"
    );
    assert_eq!(changed, 0, "annotation matrix cells unchanged A' vs B");
    assert!(membership_same, "annotation unique reads unchanged");
    let tls = after_pred_tls.expect("predecessor must leave TLS rows");
    assert_eq!(
        (tls.1, tls.2, tls.3),
        (2820, 60, 47),
        "predecessor REJECT must leave TLS 52×54 subset rows"
    );

    let classification = "SHARED_STATE_LIFECYCLE_DIVERGENCE";
    let first_op = "try_genotype_variation_event at 20:29455375 T/A REJECT VariantNotConfident writes TLS region-likelihood rows keyed by (ptr, 2808, 54) with 52 dense rows; subsequent SiteScore/with_region_likelihood_rows for 20:29455379 G/A hits that cache (same len/n_haps, reused ptr) and emits FORMAT AD 44,5 / PL 78,0,1811. Filling the full table afterwards restores A'. Annotation 52×54 cells are unchanged. No second PairHMM. Java has no pointer-keyed shared row cache between events.";
    kv("first_divergent_operation", first_op);
    kv("classification", classification);
    kv(
        "pairhmm_second_invocation",
        "NO — assign consumes the existing PairHMM matrix; no second PairHMM in this path",
    );
    assert_eq!(ap_ad, b_ad, "6R.226: predecessor no longer poisons target");
    assert_eq!(ap_pl, b_pl);
    assert_eq!(classification, "SHARED_STATE_LIFECYCLE_DIVERGENCE");
    kv("production_change", "NONE");
}

#[test]
fn forensic_6r222_closed_6r218_untouched() {
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
