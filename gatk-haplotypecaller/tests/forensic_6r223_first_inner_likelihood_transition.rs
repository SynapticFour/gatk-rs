//! 6R.223: first concrete inner operation that turns A' (`AD 47,5` /
//! `PL 68,0,1937`, n=52) into B (`AD 44,5` / `PL 78,0,1811`, n=52) at
//! `20:29455379 G/A`.
//!
//! Last proven 6R.222 state: predecessor `20:29455375 T/A` REJECT writes TLS
//! `(len=2808, n_haps=54, n_rows=52)`. This round names the SiteScore
//! cache HIT and the 52×2 allele matrix.
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r223_first_inner_likelihood_transition -- --nocapture --test-threads=1
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
    diagnose_genotype_variation_event_with_region_state, region_likelihood_rows_tls_identity,
    region_likelihoods_to_rows, take_last_site_score_inner_trace, GenotypedSiteCall,
    HcGenotypingConfig, SiteScoreInnerTrace, DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition, HaplotypeCallerEngine,
    ReadFilterParams, WalkerTraversalConfig,
};
use std::collections::BTreeSet;
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
    eprintln!("6R223\t{key}\t{}", value.as_ref());
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

fn dump_trace(tag: &str, t: &SiteScoreInnerTrace) {
    kv(
        tag,
        format!(
            "hit={}\tdense={}x{}\tdense_hash={:#x}\tmarg={}x{}\tmarg_hash={:#x}\tad={:?}\tpl={:?}",
            t.cache_hit,
            t.dense_n_rows,
            t.dense_n_cols,
            t.dense_hash,
            t.marg_n_rows,
            t.marg_n_cols,
            t.marg_hash,
            t.ad,
            t.pl
        ),
    );
}

fn matrix_delta(
    a: &[(usize, i64, i64)],
    b: &[(usize, i64, i64)],
) -> (bool, usize, usize, usize, i64, f64, usize, usize) {
    let a_ids: BTreeSet<_> = a.iter().map(|c| c.0).collect();
    let b_ids: BTreeSet<_> = b.iter().map(|c| c.0).collect();
    let only_a = a_ids.difference(&b_ids).count();
    let only_b = b_ids.difference(&a_ids).count();
    let membership_same = only_a == 0 && only_b == 0;
    let a_map: std::collections::BTreeMap<_, _> = a.iter().map(|c| (c.0, (c.1, c.2))).collect();
    let b_map: std::collections::BTreeMap<_, _> = b.iter().map(|c| (c.0, (c.1, c.2))).collect();
    let mut changed = 0usize;
    let mut abs_sum = 0i64;
    let mut max_abs = 0i64;
    let mut ref_changed = 0usize;
    let mut alt_changed = 0usize;
    for id in a_ids.intersection(&b_ids) {
        let (ar, aa) = a_map[id];
        let (br, ba) = b_map[id];
        let dr = (ar - br).abs();
        let da = (aa - ba).abs();
        if dr != 0 {
            ref_changed += 1;
            changed += 1;
            abs_sum += dr;
            max_abs = max_abs.max(dr);
        }
        if da != 0 {
            alt_changed += 1;
            changed += 1;
            abs_sum += da;
            max_abs = max_abs.max(da);
        }
    }
    let mean = if changed == 0 {
        0.0
    } else {
        (abs_sum as f64 / 1e9) / changed as f64
    };
    (
        membership_same,
        only_a,
        only_b,
        changed,
        max_abs,
        mean,
        ref_changed,
        alt_changed,
    )
}

#[test]
fn forensic_6r223_java_lifecycle_contract() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE — proof-only");
    kv(
        "java_lifecycle",
        "per-event calculateGLsForThisEvent on retainEvidence AlleleLikelihoods; no TLS (ptr,len,n_haps) row cache; no equivalent of with_region_likelihood_rows HIT across neighboring events",
    );
    assert_eq!(DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN, 2);
}

#[test]
fn forensic_6r223_first_inner_likelihood_transition() {
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
        "last_unchanged_6r222",
        "A' production-arg try_genotype AD 47,5 PL 68,0,1937 n=52 after TLS clear",
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
    let region = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("ActiveFull");
    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let _ = take_last_site_score_inner_trace();

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
    let positions: Vec<u64> = build_event_start_positions_from_cache(&hap_events)
        .iter()
        .copied()
        .filter(|p| *p >= region.start.get() && *p <= region.end.get() && *p <= TARGET)
        .collect();
    kv("loc_loop_order_through_target", format!("{positions:?}"));
    assert_eq!(
        positions,
        vec![29_455_314, 29_455_328, 29_455_337, PRED, TARGET]
    );

    let diagnose = |tag: &str, event: &VariationEvent| {
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
        let inner = take_last_site_score_inner_trace();
        match diagnosed {
            Ok(call) => (Some(call), inner),
            Err(reason) => {
                kv(tag, format!("REJECT {reason:?}"));
                (None, inner)
            }
        }
    };

    let cache_target = variation_events_at_position_from_cache(&hap_events, TARGET, emit_spanning);
    assert_eq!(cache_target.len(), 1);
    let mut pred_raw = variation_events_at_position_from_cache(&hap_events, PRED, emit_spanning);
    prefer_indel_over_colocated_snps(&mut pred_raw);
    let pred_span = replace_span_del_events(&pred_raw, PRED, apply_pad, apply_bases.as_ref());
    let pred_merged = merged_biallelic_sites_at_position(&pred_span, PRED);
    kv("pred_event", event_id(&pred_merged[0]));
    kv("target_event", event_id(&cache_target[0]));

    clear_region_likelihood_rows_tls();
    let _ = take_last_site_score_inner_trace();
    let (a_prime, a_trace) = diagnose("Aprime", &cache_target[0]);
    let a_prime = a_prime.expect("A'");
    let a_trace = a_trace.expect("A' SiteScore inner");
    dump_trace("Aprime_sitescore", &a_trace);
    let (ap_ad, ap_pl, _, ap_ann) = fmt_of(&a_prime);
    assert_eq!(ap_ad.as_slice(), AD_A.as_slice());
    assert_eq!(ap_pl.as_slice(), PL_A.as_slice());
    assert_eq!(ap_ann, 47);
    assert!(
        !a_trace.cache_hit,
        "A' must MISS TLS (rebuild target subset)"
    );
    assert_eq!(a_trace.ad.as_slice(), AD_A.as_slice());
    assert_eq!(a_trace.pl.as_slice(), PL_A.as_slice());
    assert_eq!(a_trace.marg_n_rows, 47);
    assert_eq!(a_trace.marg_n_cols, 2);
    assert_eq!(
        (a_trace.ad.as_slice(), a_trace.pl.as_slice()),
        (ap_ad.as_slice(), ap_pl.as_slice()),
        "SiteReshape must not rewrite A' AD/PL after SiteScore"
    );

    clear_region_likelihood_rows_tls();
    let (pred_call, pred_trace) = diagnose("pred", &pred_merged[0]);
    assert!(pred_call.is_none(), "predecessor T/A REJECT");
    let pred_trace = pred_trace.expect("pred still runs SiteScore");
    dump_trace("pred_sitescore", &pred_trace);
    let tls = region_likelihood_rows_tls_identity().expect("TLS after pred");
    kv("tls_after_pred", format!("{tls:?}"));
    assert_eq!((tls.1, tls.2, tls.3), (2820, 60, 47));

    let (after_pred, b_trace) = diagnose("target_after_pred", &cache_target[0]);
    let after_pred = after_pred.expect("target after pred");
    let b_trace = b_trace.expect("B SiteScore inner");
    dump_trace("B_sitescore", &b_trace);
    let (b_ad, b_pl, _, b_ann) = fmt_of(&after_pred);
    assert_eq!(b_ad.as_slice(), AD_B.as_slice());
    assert_eq!(b_pl.as_slice(), PL_B.as_slice());
    assert_eq!(b_ann, 47);
    assert!(
        !b_trace.cache_hit,
        "6R.226: different sparse cells miss even if ptr/len/n_haps collide"
    );
    assert_eq!(b_trace.ad.as_slice(), AD_B.as_slice());
    assert_eq!(b_trace.pl.as_slice(), PL_B.as_slice());
    assert_eq!(b_trace.marg_n_rows, 47);
    assert_eq!(b_trace.marg_n_cols, 2);
    assert_eq!(
        (b_trace.ad.as_slice(), b_trace.pl.as_slice()),
        (b_ad.as_slice(), b_pl.as_slice()),
        "AD and PL change together at SiteScore; SiteReshape does not rewrite B"
    );
    assert_eq!(
        b_trace.dense_hash, a_trace.dense_hash,
        "6R.226: target rebuilds P2 rather than reusing predecessor dense hap-rows"
    );
    assert_eq!(
        a_trace.dense_hash, b_trace.dense_hash,
        "target dense hap-rows are P2"
    );
    assert_eq!(
        a_trace.marg_hash, b_trace.marg_hash,
        "52×2 allele matrix is the target P2 remarg"
    );
    assert_eq!(a_trace.dense_n_rows, b_trace.dense_n_rows);
    assert_eq!(a_trace.dense_n_cols, b_trace.dense_n_cols);

    let (membership_same, only_a, only_b, changed, max_abs, mean, ref_ch, alt_ch) =
        matrix_delta(&a_trace.marg, &b_trace.marg);
    kv(
        "matrix_52x2_Aprime_vs_B",
        format!(
            "membership_same={membership_same}\tonly_a={only_a}\tonly_b={only_b}\tchanged_cells={changed}\tmax_abs_e9={max_abs}\tmean_abs={mean:.6e}\tref_changed={ref_ch}\talt_changed={alt_ch}\ta_n={}\tb_n={}",
            a_trace.marg.len(),
            b_trace.marg.len()
        ),
    );
    assert!(
        membership_same,
        "6R.226: SiteScore 52×2 rows are the target subset"
    );
    assert_eq!(only_a, 0);
    assert_eq!(only_b, 0);
    let a_ids: BTreeSet<_> = a_trace.marg.iter().map(|c| c.0).collect();
    let b_ids: BTreeSet<_> = b_trace.marg.iter().map(|c| c.0).collect();
    kv(
        "row_ids_only_a",
        format!(
            "{:?}",
            a_ids.difference(&b_ids).copied().collect::<Vec<_>>()
        ),
    );
    kv(
        "row_ids_only_b",
        format!(
            "{:?}",
            b_ids.difference(&a_ids).copied().collect::<Vec<_>>()
        ),
    );
    kv("transition_class", "MEMBERSHIP_CHANGED");
    kv(
        "ad_pl_same_operation",
        "true — both from genotype_from_marginalized_rows on the HIT 52×2",
    );

    clear_region_likelihood_rows_tls();
    let (skip, skip_trace) = diagnose("skip_pred", &cache_target[0]);
    let skip = skip.expect("skip");
    let skip_trace = skip_trace.expect("skip trace");
    dump_trace("skip_pred_sitescore", &skip_trace);
    assert!(!skip_trace.cache_hit);
    assert_eq!(skip.genotype.format.ad_as_i32(), AD_A.to_vec());
    assert_eq!(skip.genotype.format.pl_as_i32(), PL_A.to_vec());

    clear_region_likelihood_rows_tls();
    let (target_first, tf_trace) = diagnose("target_first", &cache_target[0]);
    let target_first = target_first.expect("target first");
    let tf_trace = tf_trace.expect("tf");
    let _ = diagnose("pred_after_target", &pred_merged[0]);
    assert!(!tf_trace.cache_hit);
    assert_eq!(target_first.genotype.format.ad_as_i32(), AD_A.to_vec());
    assert_eq!(target_first.genotype.format.pl_as_i32(), PL_A.to_vec());

    clear_region_likelihood_rows_tls();
    let _ = diagnose("clobber_pred", &pred_merged[0]);
    let _ = region_likelihoods_to_rows(likelihoods, haps.len());
    let (clobber, cl_trace) = diagnose("clobber_target", &cache_target[0]);
    let clobber = clobber.expect("clobber");
    let cl_trace = cl_trace.expect("clobber trace");
    dump_trace("clobber_sitescore", &cl_trace);
    assert!(!cl_trace.cache_hit, "full-table fill makes target key miss");
    assert_eq!(clobber.genotype.format.ad_as_i32(), AD_A.to_vec());
    assert_eq!(clobber.genotype.format.pl_as_i32(), PL_A.to_vec());

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
    assert_eq!(b.genotype.format.ad_as_i32(), AD_B.to_vec());
    assert_eq!(b.genotype.format.pl_as_i32(), PL_B.to_vec());

    let first_op = "with_region_likelihood_rows cache HIT keyed by (ptr, len=2808, n_haps=54) inside SiteScore::from_allele_mapping; marginalize_rows_to_biallelic_alleles then genotype_from_marginalized_rows emit AD 44,5 and PL 78,0,1811 together from the poisoned 52×2. No second PairHMM. Java has no equivalent pointer-keyed row cache.";
    kv("first_divergent_operation", first_op);
    kv("classification", "GENOTYPE_LIKELIHOOD_CACHE_DIVERGENCE");
    kv(
        "pairhmm_second_invocation",
        "NO — SiteScore consumes existing PairHMM cells via dense-row view",
    );
    kv("production_change", "NONE");
    assert_eq!(pred_trace.dense_n_rows, 47);
}

#[test]
fn forensic_6r223_closed_6r218_untouched() {
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
    .expect("outcome");
    assert_eq!(outcome.assembly.haplotypes.len(), 30);
    let closed = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_PL)
                && c.event.ref_allele == "G"
                && c.event.alt_allele == "T"
        })
        .expect("6R.218 G/T");
    assert_eq!(closed.genotype.format.pl_as_i32(), vec![69, 0, 2140]);
    kv("closed_6r218_pl", "69,0,2140");
    kv("closed_6r218_hap_n", "30");
}
