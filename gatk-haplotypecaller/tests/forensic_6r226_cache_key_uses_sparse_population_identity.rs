//! 6R.226: production TLS row cache keys exact sparse-cell identity + `n_haps`.
//! Pointer reuse of a dropped retainEvidence Vec must not return P1 for P2.
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! ONE PRODUCTION CHANGE: cache key pointer identity → sparse-cell identity.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r226_cache_key_uses_sparse_population_identity -- --nocapture --test-threads=1
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
    diagnose_genotype_variation_event_with_region_state, region_likelihood_rows_logical_identity,
    take_last_region_likelihood_rows_lookup_trace, take_last_site_score_inner_trace,
    take_region_likelihood_rows_lookup_log, with_region_likelihood_rows, HcGenotypingConfig,
    DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
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
const P2_DENSE: u64 = 0x7e3dbd79bc98a112;
const P1_DENSE: u64 = 0xfed79babc563c823;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R226\t{key}\t{}", value.as_ref());
}

fn hash_ll_rows(rows: &[gatk_haplotypecaller::genotyping::ReadLikelihoodRow]) -> u64 {
    let mut h = DefaultHasher::new();
    for row in rows {
        row.read_index.hash(&mut h);
        for v in &row.haplotype_log10_likelihoods {
            ((v * 1e9).round() as i64).hash(&mut h);
        }
    }
    h.finish()
}

#[test]
fn forensic_6r226_java_and_key_inventory() {
    kv("java_pin", JAVA_PIN);
    kv(
        "production_change",
        "ONE — with_region_likelihood_rows key is exact sparse-cell identity + n_haps",
    );
    kv(
        "key_before",
        "(as_ptr, len, n_haps) — storage identity of a temporary retainEvidence Vec",
    );
    kv(
        "key_after",
        "exact sequence of (read_index, haplotype_index, log10.to_bits()) + n_haps — not a hash",
    );
    kv(
        "java_lifecycle",
        "calculateGLsForThisEvent is per-event on retainEvidence AlleleLikelihoods; no TLS row cache",
    );
    kv(
        "collision_model",
        "exact equality of production sparse cells; a u64 hash is not used as identity",
    );
    assert_eq!(DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN, 2);
}

#[test]
fn forensic_6r226_cache_key_uses_sparse_population_identity() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }
    kv("java_pin", JAVA_PIN);
    kv("classification", "CACHE_KEY_IDENTITY_DIVERGENCE");
    kv("prior_mechanism", "CACHE_STORAGE_ALIASING (6R.224)");

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

    clear_region_likelihood_rows_tls();
    let _ = take_region_likelihood_rows_lookup_log();
    let (a_call, a_trace) = diagnose("Aprime", &cache_target[0]);
    let a_call = a_call.expect("A'");
    let a_trace = a_trace.expect("A'");
    assert!(!a_trace.cache_hit);
    assert_eq!(a_trace.dense_hash, P2_DENSE);
    kv("Aprime_dense", format!("{:#x}", a_trace.dense_hash));
    kv(
        "Aprime_format",
        format!(
            "AD={:?} PL={:?} GQ={} DP={}",
            a_call.genotype.format.ad_as_i32(),
            a_call.genotype.format.pl_as_i32(),
            a_call.genotype.format.gq.as_i32(),
            a_call.genotype.format.dp.as_i32()
        ),
    );

    let p2_cells = a_call.annotation_likelihoods.clone();
    assert_eq!(p2_cells.len(), 2820);
    let (p2_id, p2_n) = region_likelihood_rows_logical_identity(&p2_cells, haps.len());
    assert_eq!(p2_n, 60);
    assert_eq!(p2_id.len(), 2820);

    // A. Same contents, same n_haps, different allocation → same logical key → HIT.
    clear_region_likelihood_rows_tls();
    let _ = take_last_region_likelihood_rows_lookup_trace();
    let miss_hash = with_region_likelihood_rows(&p2_cells, haps.len(), |rows| hash_ll_rows(rows));
    let miss_lookup = take_last_region_likelihood_rows_lookup_trace().expect("A miss");
    assert!(!miss_lookup.hit, "first insert of P2 must MISS");
    assert_eq!(miss_hash, P2_DENSE);
    let occupy = p2_cells.clone();
    let mut distinct = Vec::with_capacity(p2_cells.len().saturating_mul(4).max(8192));
    distinct.extend_from_slice(&p2_cells);
    assert_ne!(
        distinct.as_ptr() as usize,
        p2_cells.as_ptr() as usize,
        "clone must be a different allocation"
    );
    let (clone_id, clone_n) = region_likelihood_rows_logical_identity(&distinct, haps.len());
    assert_eq!((clone_id, clone_n), (p2_id.clone(), p2_n));
    let hit_hash = with_region_likelihood_rows(&distinct, haps.len(), |rows| hash_ll_rows(rows));
    let hit_lookup = take_last_region_likelihood_rows_lookup_trace().expect("A hit");
    kv(
        "variant_A_same_contents_diff_alloc",
        format!(
            "src_ptr={:#x} clone_ptr={:#x} keys_equal=true hit={} dense={:#x}",
            p2_cells.as_ptr() as usize,
            distinct.as_ptr() as usize,
            hit_lookup.hit,
            hit_hash
        ),
    );
    assert!(
        hit_lookup.hit,
        "same sparse cells + n_haps at a different pointer must HIT"
    );
    assert_eq!(hit_hash, P2_DENSE);
    let _keep_occupy = occupy;

    // B. Different contents, same len, same n_haps, reused allocation → MISS P2.
    clear_region_likelihood_rows_tls();
    let _ = take_region_likelihood_rows_lookup_log();
    let (pred_call, pred_trace) = diagnose("pred", &pred_merged[0]);
    assert!(pred_call.is_none());
    let pred_trace = pred_trace.expect("pred");
    let pred_log = take_region_likelihood_rows_lookup_log();
    let pred_legit = pred_log.iter().filter(|e| e.legitimate_hit).count();
    let pred_alias = pred_log.iter().filter(|e| e.invalid_alias).count();
    kv(
        "pred_log",
        format!(
            "n={} legitimate_hit={} invalid_alias={} dense={:#x}",
            pred_log.len(),
            pred_legit,
            pred_alias,
            pred_trace.dense_hash
        ),
    );
    assert_eq!(pred_trace.dense_hash, P1_DENSE);
    assert!(
        pred_legit >= 1,
        "intra-event same-slice HIT must remain after the key change"
    );
    assert_eq!(pred_alias, 0);
    assert!(pred_trace.cache_hit);

    let pred_lookup = pred_trace.lookup.as_ref().expect("pred lookup");
    let (_, b_trace) = diagnose("target_after_pred", &cache_target[0]);
    let b_trace = b_trace.expect("target after pred");
    let tgt_lookup = b_trace.lookup.as_ref().expect("target lookup");
    kv(
        "variant_B_reused_alloc_diff_cells",
        format!(
            "pred_ptr={:#x} tgt_ptr={:#x} ptr_equal={} len={} n_haps={} sparse_equal={} hit={} alias={} dense={:#x}",
            pred_lookup.input_ptr,
            tgt_lookup.input_ptr,
            pred_lookup.input_ptr == tgt_lookup.input_ptr,
            tgt_lookup.input_len,
            tgt_lookup.n_haps,
            pred_lookup.input_sparse_hash == tgt_lookup.input_sparse_hash,
            tgt_lookup.hit,
            tgt_lookup.invalid_alias,
            b_trace.dense_hash
        ),
    );
    assert_eq!(pred_lookup.input_len, tgt_lookup.input_len);
    assert_eq!(pred_lookup.n_haps, tgt_lookup.n_haps);
    assert_ne!(pred_lookup.input_sparse_hash, tgt_lookup.input_sparse_hash);
    assert!(
        !tgt_lookup.hit,
        "different sparse cells must MISS even if ptr/len/n_haps collide"
    );
    assert!(!tgt_lookup.invalid_alias);
    assert_eq!(b_trace.dense_hash, P2_DENSE);
    assert_ne!(b_trace.dense_hash, P1_DENSE);

    // C. Canonical production: P2 dense rows, no alias HIT.
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
    let n_alias = log.iter().filter(|e| e.invalid_alias).count();
    let n_legit = log.iter().filter(|e| e.legitimate_hit).count();
    kv(
        "assign_region_log",
        format!(
            "n={} legitimate={} invalid_alias={}",
            log.len(),
            n_legit,
            n_alias
        ),
    );
    assert_eq!(n_alias, 0, "production assign must not alias P1 onto P2");
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
        "assign_format",
        format!(
            "AD={:?} PL={:?} GQ={} DP={}",
            site.genotype.format.ad_as_i32(),
            site.genotype.format.pl_as_i32(),
            site.genotype.format.gq.as_i32(),
            site.genotype.format.dp.as_i32()
        ),
    );
    assert_eq!(site.genotype.format.ad_as_i32(), vec![42, 5]);
    assert_eq!(site.genotype.format.pl_as_i32(), vec![84, 0, 1738]);
    assert_eq!(site.genotype.format.gq.as_i32(), 84);
    assert_eq!(site.genotype.format.dp.as_i32(), 47);

    let target_call = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based.get() == TARGET
                && c.event.ref_allele == "G"
                && c.event.alt_allele == "A"
        })
        .expect("call_region G/A");
    kv(
        "call_region_format",
        format!(
            "AD={:?} PL={:?} GQ={} DP={}",
            target_call.genotype.format.ad_as_i32(),
            target_call.genotype.format.pl_as_i32(),
            target_call.genotype.format.gq.as_i32(),
            target_call.genotype.format.dp.as_i32()
        ),
    );
    let emitted =
        try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let rec = emitted
        .iter()
        .find(|r| r.position == TARGET && r.reference == "G")
        .expect("emit G/A");
    kv(
        "vcf_format",
        format!(
            "GT={:?} AD={:?} DP={:?} GQ={:?} PL={:?} QUAL={:?}",
            rec.samples[0].gt.as_ref().map(|g| g.to_string()),
            rec.samples[0].ad.as_deref(),
            rec.samples[0].dp,
            rec.samples[0].gq,
            rec.samples[0].pl.as_deref(),
            rec.quality
        ),
    );
    assert_eq!(
        target_call.genotype.format.pl_as_i32(),
        vec![84, 0, 1738],
        "call_region must emit P2 PL, not poisoned P1 PL"
    );
}

#[test]
fn forensic_6r226_closed_6r218_untouched() {
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
