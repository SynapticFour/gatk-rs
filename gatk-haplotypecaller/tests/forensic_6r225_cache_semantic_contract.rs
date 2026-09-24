//! 6R.225: semantic contract of `with_region_likelihood_rows`.
//! Proves the cache value is an owned dense copy; the key is storage identity
//! of a short-lived subset Vec; intra-event HIT is legitimate; cross-event
//! allocator alias is invalid. Diagnostic disable / content-hash key restore P2.
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r225_cache_semantic_contract -- --nocapture --test-threads=1
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
    diagnose_genotype_variation_event_with_region_state,
    set_region_likelihood_rows_cache_diagnostic, take_last_site_score_inner_trace,
    take_region_likelihood_rows_lookup_log, HcGenotypingConfig,
    DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams,
    WalkerTraversalConfig,
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

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R225\t{key}\t{}", value.as_ref());
}

struct DiagReset;
impl Drop for DiagReset {
    fn drop(&mut self) {
        set_region_likelihood_rows_cache_diagnostic(0);
        clear_region_likelihood_rows_tls();
    }
}

#[test]
fn forensic_6r225_java_and_caller_inventory() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE — proof-only");
    kv(
        "java_lifecycle",
        "calculateGLsForThisEvent is per-event on retainEvidence AlleleLikelihoods; no TLS row cache",
    );
    kv(
        "cache_value",
        "owned Vec<ReadLikelihoodRow> { read_index: usize, read_id: String, haplotype_log10_likelihoods: Vec<f64> } — independent copied cells, no borrow of the source subset",
    );
    kv(
        "cache_key_production",
        "(as_ptr, len, n_haps) — storage identity of the current sparse slice, not sparse-cell contents",
    );
    kv(
        "production_callers",
        "with_region_likelihood_rows: SiteScore::from_allele_mapping, region_likelihoods_to_rows wrapper; wrapper callers: genotype_site_pipeline, genotype_assign, genotype_finalize, genotype_site_early_template annotation, variant_site_hc_annotations, allele_filtering, l9_dense_pileup_probe, mod.rs helpers",
    );
    kv(
        "intended_reuse_from_g1",
        "multi-allelic / multi-caller reshape of the same live sparse matrix; first call rebuilds, later calls borrow the owned dense table",
    );
    assert_eq!(DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN, 2);
}

#[test]
fn forensic_6r225_cache_semantic_contract() {
    let _reset = DiagReset;
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE — proof-only");

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
    let mut pred_raw = variation_events_at_position_from_cache(&hap_events, PRED, emit_spanning);
    prefer_indel_over_colocated_snps(&mut pred_raw);
    let pred_span = replace_span_del_events(&pred_raw, PRED, apply_pad, apply_bases.as_ref());
    let pred_merged = merged_biallelic_sites_at_position(&pred_span, PRED);

    set_region_likelihood_rows_cache_diagnostic(0);
    clear_region_likelihood_rows_tls();
    let _ = take_region_likelihood_rows_lookup_log();
    let (_, a_trace) = diagnose("Aprime", &cache_target[0]);
    let a_trace = a_trace.expect("A'");
    assert!(!a_trace.cache_hit);
    assert_eq!(a_trace.dense_hash, P2_DENSE);
    kv("Aprime_dense", format!("{:#x}", a_trace.dense_hash));

    clear_region_likelihood_rows_tls();
    let _ = take_region_likelihood_rows_lookup_log();
    let (pred_call, pred_trace) = diagnose("pred", &pred_merged[0]);
    assert!(pred_call.is_none());
    let pred_trace = pred_trace.expect("pred");
    let pred_log = take_region_likelihood_rows_lookup_log();
    let pred_legit = pred_log.iter().filter(|e| e.legitimate_hit).count();
    let pred_alias = pred_log.iter().filter(|e| e.invalid_alias).count();
    let pred_miss = pred_log.iter().filter(|e| !e.hit).count();
    kv(
        "pred_log",
        format!(
            "n={} miss={} legitimate_hit={} invalid_alias={} dense={:#x}",
            pred_log.len(),
            pred_miss,
            pred_legit,
            pred_alias,
            pred_trace.dense_hash
        ),
    );
    assert_eq!(pred_trace.dense_hash, P1_DENSE);
    assert!(
        pred_legit >= 1,
        "intra-event SiteScore HIT after region_likelihoods_to_rows on the same live subset is legitimate reuse"
    );
    assert_eq!(
        pred_alias, 0,
        "predecessor event must not alias a prior population"
    );
    assert!(pred_trace.cache_hit);

    let (_, b_trace) = diagnose("target_after_pred", &cache_target[0]);
    let b_trace = b_trace.expect("B");
    let tgt_lookup = b_trace.lookup.as_ref().expect("target lookup");
    kv(
        "target_production",
        format!(
            "hit={} legitimate={} alias={} dense={:#x} sparse={:#x}",
            tgt_lookup.hit,
            tgt_lookup.legitimate_hit,
            tgt_lookup.invalid_alias,
            b_trace.dense_hash,
            tgt_lookup.input_sparse_hash
        ),
    );
    assert!(!tgt_lookup.hit, "6R.226: P2 sparse identity misses P1 slot");
    assert!(
        !tgt_lookup.invalid_alias,
        "cross-event lookup with different sparse-cell hash must not HIT"
    );
    assert!(!tgt_lookup.legitimate_hit);
    assert_eq!(b_trace.dense_hash, P2_DENSE);
    assert_eq!(b_trace.dense_hash, a_trace.dense_hash);

    kv(
        "value_outlives_source",
        "true — cached Vec<ReadLikelihoodRow> is an owned copy; after P1 Vec drop the dense table is still memory-valid, but keyed by a recycled pointer so it is logically P1 for P2",
    );
    kv(
        "variant_C_independent_storage",
        "already the production value type; extra copying of the cached rows cannot prevent the collision because the key still aliases",
    );

    // Variant A — cache disabled.
    set_region_likelihood_rows_cache_diagnostic(1);
    clear_region_likelihood_rows_tls();
    let _ = diagnose("pred_disabled", &pred_merged[0]);
    let (_, a_dis) = diagnose("target_disabled", &cache_target[0]);
    let a_dis = a_dis.expect("disabled");
    kv(
        "variant_A_disabled",
        format!("hit={} dense={:#x}", a_dis.cache_hit, a_dis.dense_hash),
    );
    assert!(!a_dis.cache_hit);
    assert_eq!(a_dis.dense_hash, P2_DENSE, "disable restores P2");

    // Variant B — content-hash key.
    set_region_likelihood_rows_cache_diagnostic(2);
    clear_region_likelihood_rows_tls();
    let (_, pred_ck) = diagnose("pred_content_key", &pred_merged[0]);
    let pred_ck = pred_ck.expect("pred ck");
    let (_, tgt_ck) = diagnose("target_content_key", &cache_target[0]);
    let tgt_ck = tgt_ck.expect("target ck");
    let pred_key = pred_ck.lookup.as_ref().expect("pred ck lookup");
    let tgt_key = tgt_ck.lookup.as_ref().expect("tgt ck lookup");
    kv(
        "variant_B_content_key",
        format!(
            "pred_key={:#x}/{} tgt_key={:#x}/{} hit={} alias={} dense={:#x}",
            pred_key.input_ptr,
            pred_key.input_len,
            tgt_key.input_ptr,
            tgt_key.input_len,
            tgt_ck.cache_hit,
            tgt_key.invalid_alias,
            tgt_ck.dense_hash
        ),
    );
    assert_ne!(
        (pred_key.input_ptr, pred_key.input_len),
        (tgt_key.input_ptr, tgt_key.input_len),
        "content-hash keys must distinguish P1 from P2"
    );
    assert!(!tgt_key.invalid_alias);
    assert_eq!(tgt_ck.dense_hash, P2_DENSE, "content-hash key restores P2");

    set_region_likelihood_rows_cache_diagnostic(0);
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
    let log = take_region_likelihood_rows_lookup_log();
    let n_hit = log.iter().filter(|e| e.hit).count();
    let n_miss = log.iter().filter(|e| !e.hit).count();
    let n_legit = log.iter().filter(|e| e.legitimate_hit).count();
    let n_alias = log.iter().filter(|e| e.invalid_alias).count();
    kv(
        "assign_region_log",
        format!(
            "n={} miss={} hit={} legitimate={} invalid_alias={}",
            log.len(),
            n_miss,
            n_hit,
            n_legit,
            n_alias
        ),
    );
    assert!(n_miss >= 1);
    assert_eq!(
        n_alias, 0,
        "this region's production assign must not alias P1 onto P2"
    );
    kv(
        "assign_legitimate_reuse",
        format!(
            "n_legit={n_legit} — loc-loop mostly MISS (new subset Vec per event); isolated pred diagnose still has intra-event legitimate HIT (region_likelihoods_to_rows then SiteScore on the same live slice)"
        ),
    );
    let site = assigned
        .calls
        .iter()
        .find(|c| {
            c.event.start_1based.get() == TARGET
                && c.event.ref_allele == "G"
                && c.event.alt_allele == "A"
        })
        .expect("assign G/A");
    kv(
        "assign_format_consequence",
        format!(
            "AD={:?} PL={:?}",
            site.genotype.format.ad_as_i32(),
            site.genotype.format.pl_as_i32()
        ),
    );

    kv(
        "minimal_semantic_contract",
        "HIT valid iff cached dense rows were rebuilt from the same sparse-cell population as the current lookup. Pointer identity of a dropped subset Vec is not that population. Event coordinates are not required: sparse-cell contents plus n_haps determine the reshape. Intra-event repeats of the same live slice remain legitimate. Java has no equivalent cache; uncached rebuild is semantically sufficient (G1 is a performance KEEP, not a correctness KEEP).",
    );
    kv(
        "decision_case",
        "B — legitimate intra-event reuse exists; pointer identity is not semantic",
    );
    kv("classification", "CACHE_KEY_IDENTITY_DIVERGENCE");
}

#[test]
fn forensic_6r225_closed_6r218_untouched() {
    let _reset = DiagReset;
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
