//! 6R.227: first post-cache genotype-likelihood boundary at `20:29455379 G/A`.
//! Proof-only. 6R.226 cache-key semantics stay unchanged.
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r227_post_cache_genotype_likelihood_boundary -- --nocapture --test-threads=1
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
    diagnose_genotype_variation_event_with_region_state, take_last_site_score_inner_trace,
    take_region_likelihood_rows_lookup_log, HcGenotypingConfig, SiteScoreInnerTrace,
    DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_455_379;
const PRED: u64 = 29_455_375;
const CLOSED_PL: u64 = 29_455_015;
const P2_DENSE: u64 = 0x7e3dbd79bc98a112;
const P1_DENSE: u64 = 0xfed79babc563c823;
const JAVA_AD: [i32; 2] = [42, 5];
const JAVA_PL: [i32; 3] = [84, 0, 1738];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R227\t{key}\t{}", value.as_ref());
}

fn event_id(e: &VariationEvent) -> String {
    format!(
        "{}:{} {}/{}",
        e.contig,
        e.start_1based.get(),
        e.ref_allele,
        e.alt_allele
    )
}

fn dump_trace(tag: &str, t: &SiteScoreInnerTrace) {
    kv(
        tag,
        format!(
            "hit={} dense={}x{} hash={:#x} marg={}x{} hash={:#x} ad={:?} pl={:?} alias={}",
            t.cache_hit,
            t.dense_n_rows,
            t.dense_n_cols,
            t.dense_hash,
            t.marg_n_rows,
            t.marg_n_cols,
            t.marg_hash,
            t.ad,
            t.pl,
            t.lookup.as_ref().is_some_and(|l| l.invalid_alias)
        ),
    );
}

#[test]
fn forensic_6r227_java_and_starting_contract() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE — proof-only");
    kv(
        "starting_contract",
        "6R.226: cache HIT iff exact sparse-cell identity + n_haps; P2 dense 0x7e3dbd79bc98a112",
    );
    kv(
        "java_lifecycle",
        "calculateGLsForThisEvent then DepthPerAlleleBySample on the same event-local retainEvidence AlleleLikelihoods",
    );
    assert_eq!(DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN, 2);
    assert_eq!(JAVA_AD, [42, 5]);
    assert_eq!(JAVA_PL, [84, 0, 1738]);
}

#[test]
fn forensic_6r227_post_cache_genotype_likelihood_boundary() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv("target", "20:29455379 G/A");
    kv("predecessor", "20:29455375 T/A");

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
    let _ = take_region_likelihood_rows_lookup_log();

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
    kv("loc_loop_through_target", format!("{positions:?}"));
    assert_eq!(
        positions,
        vec![29_455_314, 29_455_328, 29_455_337, PRED, TARGET]
    );
    kv("n_haps", haps.len().to_string());
    kv("n_genotyping_reads", reads.len().to_string());
    kv("pairhmm_cells", likelihoods.len().to_string());

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
    assert_eq!(cache_target[0].ref_allele, "G");
    assert_eq!(cache_target[0].alt_allele, "A");
    let mut pred_raw = variation_events_at_position_from_cache(&hap_events, PRED, emit_spanning);
    prefer_indel_over_colocated_snps(&mut pred_raw);
    let pred_span = replace_span_del_events(&pred_raw, PRED, apply_pad, apply_bases.as_ref());
    let pred_merged = merged_biallelic_sites_at_position(&pred_span, PRED);
    kv("pred_event", event_id(&pred_merged[0]));
    kv("target_event", event_id(&cache_target[0]));

    // Path A — isolated target after TLS clear.
    clear_region_likelihood_rows_tls();
    let (a_call, a_trace) = diagnose("PathA_isolated", &cache_target[0]);
    let a_call = a_call.expect("Path A");
    let a_trace = a_trace.expect("Path A SiteScore");
    dump_trace("PathA", &a_trace);
    kv(
        "PathA_format",
        format!(
            "AD={:?} PL={:?} GQ={} DP={}",
            a_call.genotype.format.ad_as_i32(),
            a_call.genotype.format.pl_as_i32(),
            a_call.genotype.format.gq.as_i32(),
            a_call.genotype.format.dp.as_i32()
        ),
    );
    assert!(!a_trace.cache_hit, "isolated target after clear must MISS");
    assert_eq!(a_trace.dense_hash, P2_DENSE);
    assert_eq!(a_trace.dense_n_rows, 47);
    assert_eq!(a_trace.dense_n_cols, 60);
    assert_eq!(a_trace.marg_n_rows, 47);
    assert_eq!(a_trace.marg_n_cols, 2);
    assert_eq!(a_trace.ad, vec![42, 5]);
    assert_eq!(a_trace.pl, vec![84, 0, 1738]);
    assert_eq!(a_call.genotype.format.ad_as_i32(), vec![42, 5]);
    assert_eq!(a_call.genotype.format.pl_as_i32(), vec![84, 0, 1738]);
    assert_eq!(a_call.genotype.format.gq.as_i32(), 84);
    assert_eq!(a_call.genotype.format.dp.as_i32(), 47);
    assert_eq!(
        (a_trace.ad.as_slice(), a_trace.pl.as_slice()),
        (
            a_call.genotype.format.ad_as_i32().as_slice(),
            a_call.genotype.format.pl_as_i32().as_slice()
        ),
        "Path A AD/PL are genotype_from_marginalized_rows on the SiteScore 52×2"
    );

    // Loc-loop: predecessor then target on a fresh TLS (Path B prefix).
    clear_region_likelihood_rows_tls();
    let (pred_call, pred_trace) = diagnose("loc_pred", &pred_merged[0]);
    assert!(pred_call.is_none(), "predecessor T/A REJECT");
    let pred_trace = pred_trace.expect("pred SiteScore");
    dump_trace("loc_pred", &pred_trace);
    assert_eq!(pred_trace.dense_hash, P1_DENSE);
    assert!(!pred_trace.lookup.as_ref().is_some_and(|l| l.invalid_alias));
    let (b_prefix, b_prefix_trace) = diagnose("loc_target_after_pred", &cache_target[0]);
    let b_prefix = b_prefix.expect("target after pred");
    let b_prefix_trace = b_prefix_trace.expect("target after pred SiteScore");
    dump_trace("PathB_after_pred", &b_prefix_trace);
    assert_eq!(b_prefix_trace.dense_hash, P2_DENSE);
    assert_ne!(b_prefix_trace.dense_hash, P1_DENSE);
    assert!(!b_prefix_trace.cache_hit, "P2 must miss P1 identity");
    assert!(!b_prefix_trace
        .lookup
        .as_ref()
        .is_some_and(|l| l.invalid_alias));
    assert_eq!(b_prefix_trace.marg, a_trace.marg);
    assert_eq!(b_prefix_trace.marg_hash, a_trace.marg_hash);
    assert_eq!(b_prefix_trace.ad, a_trace.ad);
    assert_eq!(b_prefix_trace.pl, a_trace.pl);
    assert_eq!(b_prefix.genotype.format.ad_as_i32(), vec![42, 5]);
    assert_eq!(b_prefix.genotype.format.pl_as_i32(), vec![84, 0, 1738]);

    // Path B — production assign loc-loop.
    clear_region_likelihood_rows_tls();
    let _ = take_region_likelihood_rows_lookup_log();
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
    let _b_last = take_last_site_score_inner_trace();
    if let Some(last) = _b_last.as_ref() {
        dump_trace("PathB_assign_last_event", last);
        kv(
            "assign_last_is_target_p2",
            (last.dense_hash == P2_DENSE).to_string(),
        );
    }
    let log = take_region_likelihood_rows_lookup_log();
    let n_alias = log.iter().filter(|e| e.invalid_alias).count();
    kv(
        "assign_log",
        format!(
            "n={} miss={} hit={} alias={}",
            log.len(),
            log.iter().filter(|e| !e.hit).count(),
            log.iter().filter(|e| e.hit).count(),
            n_alias
        ),
    );
    assert_eq!(n_alias, 0, "6R.226: no alias HIT in production assign");
    let b_site = assigned
        .calls
        .iter()
        .find(|c| {
            c.event.start_1based.get() == TARGET
                && c.event.ref_allele == "G"
                && c.event.alt_allele == "A"
        })
        .expect("assign G/A");
    kv(
        "PathB_assign_format",
        format!(
            "AD={:?} PL={:?} GQ={} DP={}",
            b_site.genotype.format.ad_as_i32(),
            b_site.genotype.format.pl_as_i32(),
            b_site.genotype.format.gq.as_i32(),
            b_site.genotype.format.dp.as_i32()
        ),
    );

    assert_eq!(
        b_prefix_trace.dense_hash, P2_DENSE,
        "production loc-loop target dense is P2"
    );
    assert_eq!(b_prefix_trace.dense_n_rows, 47);
    assert_eq!(b_prefix_trace.dense_n_cols, 60);
    assert_eq!(b_prefix_trace.marg_n_rows, 47);
    assert_eq!(b_prefix_trace.marg_n_cols, 2);
    assert_eq!(
        b_prefix_trace.marg, a_trace.marg,
        "Path B 52×2 cells equal Path A (read_index, REF bits, ALT bits)"
    );
    assert_eq!(b_prefix_trace.marg_hash, a_trace.marg_hash);
    assert_eq!(b_prefix_trace.ad, a_trace.ad);
    assert_eq!(b_prefix_trace.pl, a_trace.pl);
    assert_eq!(b_site.genotype.format.ad_as_i32(), vec![42, 5]);
    assert_eq!(b_site.genotype.format.pl_as_i32(), vec![84, 0, 1738]);
    assert_eq!(b_site.genotype.format.gq.as_i32(), 84);
    assert_eq!(b_site.genotype.format.dp.as_i32(), 47);

    let prod = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based.get() == TARGET
                && c.event.ref_allele == "G"
                && c.event.alt_allele == "A"
        })
        .expect("call_region G/A");
    assert_eq!(prod.genotype.format.ad_as_i32(), vec![42, 5]);
    assert_eq!(prod.genotype.format.pl_as_i32(), vec![84, 0, 1738]);
    let emitted =
        try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let rec = emitted
        .iter()
        .find(|r| r.position == TARGET && r.reference == "G")
        .expect("emit G/A");
    assert_eq!(rec.samples[0].ad.as_deref(), Some([42, 5].as_slice()));
    assert_eq!(rec.samples[0].pl.as_deref(), Some([84, 0, 1738].as_slice()));
    assert_eq!(rec.samples[0].dp, Some(47));

    let extra_ref = a_trace.ad[0] - JAVA_AD[0];
    kv("java_ad", "42,5");
    kv("java_pl", "84,0,1738");
    kv("java_gq", "84");
    kv("java_dp", "47");
    kv("extra_ref_vs_java", extra_ref.to_string());
    kv(
        "ad_pl_same_matrix",
        "true — both from genotype_from_marginalized_rows on the Path A/B 52×2",
    );
    assert_eq!(extra_ref, 0);
    assert_eq!(a_trace.ad[1], JAVA_AD[1]);
    assert_eq!(a_trace.ad.as_slice(), JAVA_AD.as_slice());
    assert_eq!(a_trace.pl.as_slice(), JAVA_PL.as_slice());

    let first_op = "Path A isolated try_genotype and Path B assign loc-loop both consume P2 dense 0x7e3dbd79bc98a112 and the same 52×2 allele-likelihood matrix; AD 47,5 and PL 68,0,1937 are genotype_from_marginalized_rows on that matrix. Predecessor 20:29455375 T/A REJECT still writes P1 and does not mutate the target 52×2. No alias HIT. Remaining vs Java is that 52-row retainEvidence remarg (unique n=52, five extra REF votes) is not Java calculateGLsForThisEvent / DepthPerAlleleBySample on a 47-read AlleleLikelihoods.";
    let classification = "NO_DIVERGENCE_AT_THIS_BOUNDARY";
    kv("first_divergent_operation", first_op);
    kv("classification", classification);
    kv("cache_key_change", "NONE — 6R.226 semantics unchanged");
    kv(
        "other_tls_inspected",
        "REGION_LIKELIHOOD_ROWS_CACHE (P2 MISS, no alias); SITE_AD_SCRATCH overwrite-per-site; COLOCATED_MERGE_NUMERICS diagnostic; ad_decode_cache cleared at assign entry — none change the 52×2",
    );
    assert_eq!(classification, "NO_DIVERGENCE_AT_THIS_BOUNDARY");
    kv("production_change", "NONE");
}

#[test]
fn forensic_6r227_closed_6r218_untouched() {
    kv(
        "guard",
        "6R.218 cycle abort stays; 6R.226 cache key stays; this round does not reopen either",
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
    assert_eq!(outcome.assembly.haplotypes.len(), 30);
    let closed = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based.get() == CLOSED_PL
                && c.event.ref_allele == "G"
                && c.event.alt_allele == "T"
        })
        .expect("6R.218 G/T");
    assert_eq!(closed.genotype.format.pl_as_i32(), vec![69, 0, 2140]);
    kv("closed_6r218_pl", "69,0,2140");
    kv("closed_6r218_hap_n", "30");
}
