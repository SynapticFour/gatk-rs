//! 6R.224: why predecessor `20:29455375 T/A` and target `20:29455379 G/A`
//! produce the same TLS row-cache key while their logical populations differ.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE. Cache key / invalidation / PairHMM unchanged.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r224_tls_cache_key_identity -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{
    build_event_start_positions_from_cache, build_per_haplotype_variation_events,
    merged_biallelic_sites_at_position, prefer_indel_over_colocated_snps,
    variation_events_at_position_from_cache, VariationEvent,
};
use gatk_haplotypecaller::genotyping::ReadLikelihoodRow;
use gatk_haplotypecaller::hc_allele_mapping::replace_span_del_events;
use gatk_haplotypecaller::hc_genotyping_engine::{
    clear_region_likelihood_rows_tls, diagnose_genotype_variation_event_with_region_state,
    region_likelihood_rows_tls_identity, region_likelihoods_to_rows_uncached_pub,
    take_last_region_likelihood_rows_lookup_trace, take_last_site_score_inner_trace,
    with_region_likelihood_rows, HcGenotypingConfig, RegionLikelihoodRowsLookupTrace,
    DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams,
    WalkerTraversalConfig,
};
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_455_379;
const PRED: u64 = 29_455_375;
const CLOSED_PL: u64 = 29_455_015;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R224\t{key}\t{}", value.as_ref());
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

fn hash_ll_rows(rows: &[ReadLikelihoodRow]) -> u64 {
    let mut h = DefaultHasher::new();
    for row in rows {
        row.read_index.hash(&mut h);
        for v in &row.haplotype_log10_likelihoods {
            ((v * 1e9).round() as i64).hash(&mut h);
        }
    }
    h.finish()
}

fn dump_lookup(tag: &str, t: &RegionLikelihoodRowsLookupTrace) {
    kv(
        tag,
        format!(
            "hit={}\tptr={:#x}\tlen={}\tn_haps={}\tsparse={:#x}\tread_set={:#x}\tn_unique={}\tstored_ptr={:?}\tstored_len={:?}\tstored_n_haps={:?}\tstored_n_rows={:?}",
            t.hit,
            t.input_ptr,
            t.input_len,
            t.n_haps,
            t.input_sparse_hash,
            t.input_read_set_hash,
            t.input_n_unique_reads,
            t.stored_ptr.map(|p| format!("{p:#x}")),
            t.stored_len,
            t.stored_n_haps,
            t.stored_n_rows
        ),
    );
}

#[test]
fn forensic_6r224_java_cache_contract() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE — proof-only");
    kv(
        "java_lifecycle",
        "calculateGLsForThisEvent is per-event on that event's retainEvidence AlleleLikelihoods; Java has no TLS (ptr,len,n_haps) dense-row cache and does not memoize likelihood rows across neighboring events",
    );
    kv(
        "rust_cache_contract",
        "with_region_likelihood_rows TLS keyed by (as_ptr, len, n_haps); intended for multi-allelic repeats of the same live sparse slice; production never clears the slot",
    );
    assert_eq!(DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN, 2);
}

#[test]
fn forensic_6r224_tls_cache_key_identity() {
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
    let _ = take_last_region_likelihood_rows_lookup_trace();

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
    kv("pairhmm_n_cells", likelihoods.len().to_string());
    kv("n_haps", haps.len().to_string());
    assert_eq!(haps.len(), 60);

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
        let lookup = inner.as_ref().and_then(|t| t.lookup.clone());
        match diagnosed {
            Ok(call) => (Some(call), inner, lookup),
            Err(reason) => {
                kv(tag, format!("REJECT {reason:?}"));
                (None, inner, lookup)
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
    assert_ne!(
        event_id(&pred_merged[0]),
        event_id(&cache_target[0]),
        "predecessor and target are distinct events"
    );

    clear_region_likelihood_rows_tls();
    let _ = take_last_site_score_inner_trace();
    let _ = take_last_region_likelihood_rows_lookup_trace();
    let (a_prime, a_trace, a_lookup) = diagnose("Aprime", &cache_target[0]);
    let a_prime = a_prime.expect("A'");
    let a_trace = a_trace.expect("A' SiteScore inner");
    let a_lookup = a_lookup.expect("A' lookup");
    dump_lookup("Aprime_lookup", &a_lookup);
    kv(
        "Aprime_dense",
        format!(
            "{}x{} hash={:#x}",
            a_trace.dense_n_rows, a_trace.dense_n_cols, a_trace.dense_hash
        ),
    );
    assert!(
        !a_lookup.hit,
        "isolated target after TLS clear must MISS and rebuild P2"
    );
    assert_eq!(a_lookup.input_len, 2820);
    assert_eq!(a_lookup.n_haps, 60);
    assert_eq!(a_lookup.input_n_unique_reads, 47);
    let p2_uncached =
        region_likelihoods_to_rows_uncached_pub(&a_prime.annotation_likelihoods, haps.len());
    kv(
        "P2_annotation_uncached",
        format!(
            "n_cells={} n_rows={} hash={:#x}",
            a_prime.annotation_likelihoods.len(),
            p2_uncached.len(),
            hash_ll_rows(&p2_uncached)
        ),
    );
    assert_eq!(
        hash_ll_rows(&p2_uncached),
        a_trace.dense_hash,
        "annotation cells are the target P2 hap-rows"
    );

    clear_region_likelihood_rows_tls();
    let (pred_call, pred_trace, pred_lookup) = diagnose("pred", &pred_merged[0]);
    assert!(pred_call.is_none(), "predecessor T/A REJECT");
    let pred_trace = pred_trace.expect("pred SiteScore");
    let pred_lookup = pred_lookup.expect("pred lookup");
    dump_lookup("pred_lookup", &pred_lookup);
    kv(
        "pred_dense",
        format!(
            "{}x{} hash={:#x}",
            pred_trace.dense_n_rows, pred_trace.dense_n_cols, pred_trace.dense_hash
        ),
    );
    let tls = region_likelihood_rows_tls_identity().expect("TLS after pred");
    kv("tls_after_pred", format!("{tls:?}"));
    assert_eq!((tls.1, tls.2, tls.3), (2820, 60, 47));
    assert_eq!(pred_lookup.input_len, 2820);
    assert_eq!(pred_lookup.n_haps, 60);
    assert_eq!(pred_lookup.input_n_unique_reads, 47);
    assert_eq!(pred_trace.dense_hash, 0xfed79babc563c823);
    kv(
        "cache_lifetime",
        "thread_local RefCell; production never calls clear_region_likelihood_rows_tls; slot survives events/regions until a miss overwrites it",
    );

    let (after_pred, b_trace, tgt_lookup) = diagnose("target_after_pred", &cache_target[0]);
    let after_pred = after_pred.expect("target after pred");
    let b_trace = b_trace.expect("B SiteScore");
    let tgt_lookup = tgt_lookup.expect("target lookup");
    dump_lookup("target_lookup", &tgt_lookup);
    kv(
        "B_dense",
        format!(
            "{}x{} hash={:#x}",
            b_trace.dense_n_rows, b_trace.dense_n_cols, b_trace.dense_hash
        ),
    );

    assert!(
        !tgt_lookup.hit,
        "6R.226: different sparse cells must miss even if ptr/len/n_haps collide"
    );
    assert_eq!(
        (tgt_lookup.input_len, tgt_lookup.n_haps),
        (pred_lookup.input_len, pred_lookup.n_haps)
    );
    assert_eq!(tgt_lookup.stored_len, Some(2820));
    assert_eq!(tgt_lookup.stored_n_haps, Some(60));
    assert_ne!(
        tgt_lookup.input_sparse_hash, pred_lookup.input_sparse_hash,
        "same key, different sparse-cell contents"
    );
    assert_ne!(
        tgt_lookup.input_read_set_hash, pred_lookup.input_read_set_hash,
        "same key, different read-index membership"
    );
    assert_eq!(
        tgt_lookup.input_sparse_hash, a_lookup.input_sparse_hash,
        "target lookup input is P2, not P1"
    );
    assert_eq!(
        tgt_lookup.input_read_set_hash, a_lookup.input_read_set_hash,
        "target lookup read-set is P2"
    );
    assert_eq!(
        b_trace.dense_hash, a_trace.dense_hash,
        "MISS rebuilds target P2"
    );
    assert_ne!(
        b_trace.dense_hash, pred_trace.dense_hash,
        "returned population is not predecessor P1"
    );
    assert_eq!(a_trace.dense_hash, 0x7e3dbd79bc98a112);

    kv(
        "missing_key_dimension",
        "sparse-cell contents + read-index membership + event identity; key uses the backing pointer of a dropped per-event subset Vec (filter_likelihoods_for_variant .cloned().collect()), which the allocator reuses for the next len=2808 subset",
    );
    kv(
        "ptr_points_to",
        "data pointer of the owned retainEvidence subset Vec (Cow::Owned / into_owned), not the region PairHMM table, not a persistent scratch buffer",
    );
    kv(
        "len_n_haps_protection",
        "none here: both events are 52 unique reads × 54 haplotypes = 2808 sparse cells",
    );
    kv("classification", "CACHE_STORAGE_ALIASING");

    // Experiment A: after a fresh predecessor insert (TLS still P1 at 2808),
    // a distinct allocation of P2 contents must miss and rebuild P2.
    clear_region_likelihood_rows_tls();
    let _ = diagnose("pred_for_expA", &pred_merged[0]);
    let tls_p1 = region_likelihood_rows_tls_identity().expect("P1 TLS");
    assert_eq!((tls_p1.1, tls_p1.2, tls_p1.3), (2820, 60, 47));
    let p2_cells = a_prime.annotation_likelihoods.clone();
    let occupy_recycled_slot = p2_cells.clone();
    let mut distinct = Vec::with_capacity(p2_cells.len().saturating_mul(4).max(8192));
    distinct.extend_from_slice(&p2_cells);
    kv(
        "expA_ptrs",
        format!(
            "p1_tls={:#x} occupy={:#x} distinct={:#x}",
            tls_p1.0,
            occupy_recycled_slot.as_ptr() as usize,
            distinct.as_ptr() as usize
        ),
    );
    assert_ne!(
        distinct.as_ptr() as usize,
        p2_cells.as_ptr() as usize,
        "diagnostic copy must be a different allocation"
    );
    let _keep_occupy = occupy_recycled_slot;
    let _ = take_last_region_likelihood_rows_lookup_trace();
    let rebuilt = with_region_likelihood_rows(&distinct, haps.len(), |rows| {
        (rows.len(), hash_ll_rows(rows))
    });
    let exp_a = take_last_region_likelihood_rows_lookup_trace().expect("exp A lookup");
    dump_lookup("expA_distinct_alloc", &exp_a);
    assert!(!exp_a.hit, "distinct ptr must miss the P1 key");
    assert_eq!(exp_a.stored_len, Some(2820));
    assert_eq!(exp_a.stored_n_haps, Some(60));
    assert_eq!(rebuilt.0, 47);
    assert_eq!(rebuilt.1, a_trace.dense_hash, "MISS rebuilds P2");
    kv("experiment_A_distinct_ptr_restores_P2", "true");

    // Experiment B is the production loc-loop itself: same ptr, replaced logical contents, HIT returns P1.
    kv(
        "experiment_B_same_ptr_different_contents",
        "true — target lookup input_ptr equals pred ptr while sparse/read-set hashes are P2",
    );

    let _ = after_pred;
}

#[test]
fn forensic_6r224_closed_6r218_untouched() {
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
