//! 6R.234 live: Java/Rust share FNV `c7acc50dfb9f9ecc` but not k-best
//! path state (18 vs 19 edges). The extra Rust sink-split edge is
//! log10(63/78) and moves the carrier from Java rank 116 / −2.64786702
//! to Rust rank 129 / −2.74062107. K=128 is a downstream cutoff.
//! Do not raise K. 6R.226–6R.233 unchanged. PRODUCTION CHANGE: NONE.
//! Skipped unless `HOLDOUT_6R234=1`.
//!
//! ```text
//! HOLDOUT_6R234=1 cargo test -p gatk-haplotypecaller --test holdout_6r234_c_carrier_kbest_rank -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{
    build_per_haplotype_variation_events, merged_biallelic_sites_at_position,
    prefer_indel_over_colocated_snps, variation_events_at_position_from_cache,
};
use gatk_haplotypecaller::hc_allele_mapping::replace_span_del_events;
use gatk_haplotypecaller::hc_genotyping_engine::{
    assign_genotype_likelihoods_for_region, clear_region_likelihood_rows_tls,
    diagnose_genotype_variation_event, diagnose_genotype_variation_event_with_region_state,
    java_alignment_read_overlaps_interval, region_likelihood_rows_tls_identity,
    region_likelihoods_to_rows, set_region_likelihood_rows_cache_diagnostic,
    take_last_site_score_inner_trace, HcGenotypingConfig,
};
use gatk_haplotypecaller::variant_site_hc_annotations::coverage_evidence_count;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const CHR20_INTERVAL: &str = "20:29455000-29456500";
const DP_INTERVAL: &str = "2:92305500-92305850";
const TG_INTERVAL: &str = "2:92307200-92307550";
const MID_INTERVAL: &str = "2:92316200-92316580";
const MID_B_INTERVAL: &str = "2:92317000-92319000";
const P12_POST_INTERVAL: &str = "2:92318150-92319220";
const HET_TAIL_INTERVAL: &str = "2:92324900-92325400";
const CHR20_BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const P12_BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_455_379;
const PRED: u64 = 29_455_375;
const CLOSED_218: u64 = 29_455_015;
const SIB: u64 = 29_455_019;
const CLOSED_DP: u64 = 92_305_634;
const CLOSED_AG: u64 = 92_305_635;
const CLOSED_AC: u64 = 92_305_716;
const CLOSED_TTC: u64 = 92_307_324;
const CLOSED_TG: u64 = 92_307_333;
const CLOSED_CT: u64 = 92_307_359;
const CLOSED_CA: u64 = 92_307_403;
const CLOSED_DEL: u64 = 92_316_347;
const CLOSED_HOM_ALT: u64 = 92_316_296;
const CLOSED_ONE_READ: u64 = 92_316_416;
const CLOSED_MID_B: u64 = 92_317_399;
const CLOSED_POST: u64 = 92_318_199;
const CLOSED_HET: u64 = 92_325_193;
const CLOSED_SIB: u64 = 92_325_205;
const CLOSED_WEAK: u64 = 92_325_268;
const MARGIN: i32 = 2;
const JAVA_POST_TSV: &str = include_str!("forensic_6r228_java_post_retain.tsv");
const RUST_ONLY_PIN: &[(&str, u16)] = &[];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R234\t{key}\t{}", value.as_ref());
}

struct DiagReset;
impl Drop for DiagReset {
    fn drop(&mut self) {
        set_region_likelihood_rows_cache_diagnostic(0);
        clear_region_likelihood_rows_tls();
    }
}

fn unique_indices(likelihoods: &[gatk_haplotypecaller::RegionReadLikelihood]) -> BTreeSet<usize> {
    likelihoods.iter().map(|c| c.read_index.get()).collect()
}

fn info_i32(info: &[InfoValue], key: &str) -> Option<i32> {
    for v in info {
        if let InfoValue::Integer(k, xs) = v {
            if k == key {
                return xs.first().copied();
            }
        }
    }
    None
}

#[test]
fn holdout_6r234_c_carrier_kbest_rank() {
    if std::env::var("HOLDOUT_6R234").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R234=1");
        return;
    }
    let _reset = DiagReset;
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let chr20_bam = root.join(CHR20_BAM_REL);
    let p12_bam = root.join(P12_BAM_REL);
    if !ref_fasta.is_file() || !chr20_bam.is_file() || !p12_bam.is_file() {
        eprintln!("skip: missing BAM/REF");
        return;
    }
    kv(
        "java_pin",
        "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77",
    );
    kv(
        "production_change",
        "NONE — proof-only: same FNV, Java 18-edge vs Rust 19-edge k-best state; extra sink edge log10(63/78) explains rank 116→129. Do not raise K",
    );
    kv("target", "20:29455379 G/A");
    kv(
        "predecessor",
        "6R.231 READ_INTERVAL_INPUT_DIVERGENCE CLOSED",
    );
    kv("classification", "HAPLOTYPE_CONSTRUCTION_DIVERGENCE");

    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let walk = |interval: &str, bam: &Path| {
        let specs = parse_intervals_cli_string(&dict, interval).expect("interval");
        traverse_assembly_region_walker(
            &dict,
            &specs,
            &ref_fasta,
            bam,
            &ReadFilterParams::gatk_standard_hc(),
            &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
        )
        .expect("walk")
    };
    let chr20_regions = flatten_assembly_regions(&walk(CHR20_INTERVAL, &chr20_bam));
    let covering_target = chr20_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("chr20 target");
    let covering_218 = chr20_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_218
                && r.end.get() >= CLOSED_218
        })
        .expect("6R.218 ActiveFull");
    let target_outcome = HaplotypeCallerEngine::call_region(
        covering_target,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("outcome");
    let target_call = target_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(TARGET)
                && c.event.ref_allele == "G"
                && c.event.alt_allele == "A"
        })
        .expect("G/A");
    assert!(!target_call.post_merge_unused_alt_subset);
    assert_eq!(target_call.genotype.format.ad_as_i32(), vec![42, 5]);
    assert_eq!(target_call.genotype.format.pl_as_i32(), vec![84, 0, 1738]);
    assert_eq!(target_call.genotype.format.dp.as_i32(), 47);
    let emitted = try_emit_call_region_variants(
        covering_target,
        &target_outcome,
        "SAMPLE",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("emit");
    let target_rec = emitted
        .into_iter()
        .find(|r| r.position == TARGET && r.reference == "G")
        .expect("emit G/A");
    assert_eq!(
        target_rec.samples[0].ad.as_deref(),
        Some([42, 5].as_slice())
    );
    assert_eq!(target_rec.samples[0].dp, Some(47));
    assert_eq!(info_i32(&target_rec.info, "DP"), Some(47));
    kv("java_ad", "42,5");
    kv("rust_ad", "42,5");
    kv("matched_format", "GT=0/1 AD=42,5 DP=47 matches Java");
    kv(
        "unmatched",
        "same FNV; Java 18 edges score -2.64786702 rank 116; Rust 19 edges score -2.74062107 rank 129; prefix 0..16 identical",
    );
    kv(
        "first_divergent_arrow",
        "ASSEMBLY_KBEST_STATE_IDENTITY_DIVERGENCE: same FNV, different edge walk (18 vs 19); extra Rust sink-split penalty -0.09275. K=128 is downstream. Do not raise K. 6R.226–6R.233 unchanged",
    );

    let java_ids: BTreeSet<(String, u16)> = JAVA_POST_TSV
        .lines()
        .skip(1)
        .filter(|l| !l.is_empty())
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f[0].to_string(), f[1].parse::<u16>().unwrap())
        })
        .collect();
    assert_eq!(java_ids.len(), 47);
    let pairhmm_idx: BTreeSet<usize> = target_outcome
        .read_likelihoods
        .iter()
        .map(|c| c.read_index.get())
        .collect();
    let rust_ids: BTreeSet<(String, u16)> = target_outcome
        .genotyping_reads
        .iter()
        .enumerate()
        .filter(|(i, r)| {
            pairhmm_idx.contains(i)
                && java_alignment_read_overlaps_interval(r, TARGET, TARGET, MARGIN)
        })
        .map(|(_, r)| (String::from_utf8_lossy(r.qname()).into_owned(), r.flags()))
        .collect();
    let rust_only: BTreeSet<_> = rust_ids.difference(&java_ids).cloned().collect();
    let java_only: BTreeSet<_> = java_ids.difference(&rust_ids).cloned().collect();
    assert_eq!(rust_ids.len(), 47);
    assert_eq!(rust_ids.intersection(&java_ids).count(), 47);
    assert_eq!(rust_only.len(), 0);
    assert_eq!(java_only.len(), 0);
    let pin: BTreeSet<(String, u16)> = RUST_ONLY_PIN
        .iter()
        .map(|(q, f)| ((*q).to_string(), *f))
        .collect();
    assert_eq!(rust_only, pin);
    kv("java_membership", "47");
    kv("rust_membership", "47");
    kv("intersection", "47");
    kv("rust_only", "0");
    kv("java_only", "0");

    let (full_ref, full_pad) = target_outcome.assembly.event_map_reference();
    let apply_bases = target_outcome.assembly.apply_bases_shared();
    let apply_pad = target_outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.map(|g| g.start_1based()))
        .unwrap_or(full_pad);
    let hap_events = build_per_haplotype_variation_events(
        &target_outcome.assembly.haplotypes,
        full_ref,
        full_pad,
        target_outcome.assembly.max_mnp_distance(),
        covering_target.contig.as_str(),
    );
    let cache_only = variation_events_at_position_from_cache(
        &hap_events,
        TARGET,
        !HcGenotypingConfig::strict_java().disable_spanning_event_genotyping,
    );
    assert_eq!(cache_only.len(), 1);
    let diagnosed = diagnose_genotype_variation_event(
        &cache_only[0],
        &target_outcome.read_likelihoods,
        &target_outcome.genotyping_reads,
        &covering_target.reads,
        Some(covering_target.reads.as_slice()),
        &target_outcome.assembly.haplotypes,
        apply_bases.as_ref(),
        apply_pad,
        full_ref,
        full_pad,
        covering_target.start.get(),
        covering_target.end.get(),
        target_outcome.assembly.max_mnp_distance(),
        &HcGenotypingConfig::strict_java(),
    )
    .expect("diagnose");
    let diagnosed = diagnosed.expect("try_genotype Ok");
    let path_a = take_last_site_score_inner_trace().expect("Path A SiteScore");
    assert_eq!(path_a.dense_hash, 0x7e3dbd79bc98a112);
    assert_eq!(path_a.ad, vec![42, 5]);
    assert_eq!(path_a.pl, vec![84, 0, 1738]);
    assert_eq!(path_a.marg_n_rows, 47);
    assert_eq!(path_a.marg_n_cols, 2);
    assert_eq!(diagnosed.genotype.format.ad_as_i32(), vec![42, 5]);
    assert_eq!(diagnosed.genotype.format.pl_as_i32(), vec![84, 0, 1738]);
    assert_eq!(
        diagnosed.genotype.format.ad_as_i32(),
        target_call.genotype.format.ad_as_i32(),
        "6R.226: isolated try_genotype matches production FORMAT P2 remarg"
    );
    kv("diagnose_try_genotype_ad", "42,5");
    kv("diagnose_try_genotype_pl", "84,0,1738");
    let stored = target_outcome.assembly.variation_events();
    let with_state = diagnose_genotype_variation_event_with_region_state(
        &cache_only[0],
        &target_outcome.read_likelihoods,
        &target_outcome.genotyping_reads,
        &covering_target.reads,
        Some(covering_target.reads.as_slice()),
        &target_outcome.assembly.haplotypes,
        apply_bases.as_ref(),
        apply_pad,
        full_ref,
        full_pad,
        covering_target.start.get(),
        covering_target.end.get(),
        target_outcome.assembly.max_mnp_distance(),
        &HcGenotypingConfig::strict_java(),
        stored,
        Some(&hap_events),
    )
    .expect("diagnose with state");
    let with_state = with_state.expect("try_genotype with state Ok");
    assert_eq!(with_state.genotype.format.ad_as_i32(), vec![42, 5]);
    assert_eq!(with_state.genotype.format.pl_as_i32(), vec![84, 0, 1738]);
    let assigned = assign_genotype_likelihoods_for_region(
        &target_outcome.read_likelihoods,
        &target_outcome.genotyping_reads,
        &covering_target.reads,
        Some(covering_target.reads.as_slice()),
        &target_outcome.assembly.haplotypes,
        apply_bases.as_ref(),
        apply_pad,
        full_ref,
        full_pad,
        covering_target.start.get(),
        covering_target.end.get(),
        covering_target.contig.as_str(),
        target_outcome.assembly.max_mnp_distance(),
        &HcGenotypingConfig::strict_java(),
        stored,
        &[],
    )
    .expect("assign");
    let assigned_site = assigned
        .calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(TARGET)
                && c.event.ref_allele == "G"
                && c.event.alt_allele == "A"
        })
        .expect("assign G/A");
    assert_eq!(assigned_site.genotype.format.ad_as_i32(), vec![42, 5]);
    assert_eq!(assigned_site.genotype.format.pl_as_i32(), vec![84, 0, 1738]);
    kv("assign_ad", "42,5");
    kv("assign_pl", "84,0,1738");

    clear_region_likelihood_rows_tls();
    let emit_spanning = !HcGenotypingConfig::strict_java().disable_spanning_event_genotyping;
    let cache_target = variation_events_at_position_from_cache(&hap_events, TARGET, emit_spanning);
    let mut pred_raw = variation_events_at_position_from_cache(&hap_events, PRED, emit_spanning);
    prefer_indel_over_colocated_snps(&mut pred_raw);
    let pred_span = replace_span_del_events(&pred_raw, PRED, apply_pad, apply_bases.as_ref());
    let pred_merged = merged_biallelic_sites_at_position(&pred_span, PRED);
    assert_eq!(pred_merged.len(), 1);
    assert_eq!(pred_merged[0].ref_allele, "T");
    assert_eq!(pred_merged[0].alt_allele, "A");
    let diagnose = |event: &gatk_haplotypecaller::event_map::VariationEvent| {
        diagnose_genotype_variation_event_with_region_state(
            event,
            &target_outcome.read_likelihoods,
            &target_outcome.genotyping_reads,
            &covering_target.reads,
            Some(covering_target.reads.as_slice()),
            &target_outcome.assembly.haplotypes,
            apply_bases.as_ref(),
            apply_pad,
            full_ref,
            full_pad,
            covering_target.start.get(),
            covering_target.end.get(),
            target_outcome.assembly.max_mnp_distance(),
            &HcGenotypingConfig::strict_java(),
            stored,
            Some(&hap_events),
        )
        .expect("diagnose")
    };
    let pred_res = diagnose(&pred_merged[0]);
    assert!(pred_res.is_err(), "29455375 T/A must REJECT");
    let pred_inner = take_last_site_score_inner_trace().expect("pred inner");
    let pred_lookup = pred_inner.lookup.as_ref().expect("pred SiteScore lookup");
    let tls = region_likelihood_rows_tls_identity().expect("TLS after pred");
    assert_eq!((tls.1, tls.2, tls.3), (2820, 60, 47));
    assert_eq!(pred_lookup.input_len, 2820);
    assert_eq!(pred_lookup.n_haps, 60);
    let after_pred = diagnose(&cache_target[0]).expect("target after pred");
    let b_inner = take_last_site_score_inner_trace().expect("B inner");
    let tgt_lookup = b_inner.lookup.as_ref().expect("target SiteScore lookup");
    assert!(
        !b_inner.cache_hit,
        "6R.226: sparse identity misses across events"
    );
    assert!(!tgt_lookup.hit);
    assert_eq!(
        (tgt_lookup.input_len, tgt_lookup.n_haps),
        (pred_lookup.input_len, pred_lookup.n_haps)
    );
    assert_ne!(tgt_lookup.input_sparse_hash, pred_lookup.input_sparse_hash);
    assert_ne!(
        tgt_lookup.input_read_set_hash,
        pred_lookup.input_read_set_hash
    );
    assert_ne!(b_inner.dense_hash, pred_inner.dense_hash);
    assert_eq!(b_inner.ad, vec![42, 5]);
    assert_eq!(b_inner.pl, vec![84, 0, 1738]);
    assert_eq!(b_inner.marg_n_rows, 47);
    assert_eq!(b_inner.marg_n_cols, 2);
    assert_eq!(after_pred.genotype.format.ad_as_i32(), vec![42, 5]);
    assert_eq!(after_pred.genotype.format.pl_as_i32(), vec![84, 0, 1738]);
    kv("classification_6r226_cache", "CLOSED — P2 dense, no alias");
    kv("classification", "HAPLOTYPE_CONSTRUCTION_DIVERGENCE");
    assert_eq!(
        b_inner.marg, path_a.marg,
        "Path B after predecessor has Path A 52×2 cells"
    );
    clear_region_likelihood_rows_tls();
    let skip = diagnose(&cache_target[0]).expect("skip pred");
    let skip_inner = take_last_site_score_inner_trace().expect("skip inner");
    assert!(!skip_inner.cache_hit);
    assert_eq!(skip.genotype.format.ad_as_i32(), vec![42, 5]);
    assert_eq!(skip.genotype.format.pl_as_i32(), vec![84, 0, 1738]);
    let _ = diagnose(&pred_merged[0]);
    let _ = region_likelihoods_to_rows(
        &target_outcome.read_likelihoods,
        target_outcome.assembly.haplotypes.len(),
    );
    let clobber = diagnose(&cache_target[0]).expect("clobber");
    assert_eq!(clobber.genotype.format.ad_as_i32(), vec![42, 5]);
    assert_eq!(clobber.genotype.format.pl_as_i32(), vec![84, 0, 1738]);
    kv("pred_reject_tls", "2820,60,47");
    kv("skip_pred_stays_Aprime", "true");
    kv("clobber_restores_Aprime", "true");
    kv("sitescore_hit_is_first_op", "true");
    kv(
        "key_aliasing",
        "ptr+len+n_haps equal; sparse/read-set unequal",
    );

    set_region_likelihood_rows_cache_diagnostic(1);
    clear_region_likelihood_rows_tls();
    let _ = diagnose(&pred_merged[0]);
    let disabled = diagnose(&cache_target[0]).expect("disabled");
    let dis_inner = take_last_site_score_inner_trace().expect("disabled inner");
    assert!(!dis_inner.cache_hit);
    assert_eq!(disabled.genotype.format.ad_as_i32(), vec![42, 5]);
    assert_eq!(disabled.genotype.format.pl_as_i32(), vec![84, 0, 1738]);
    kv("variant_A_disabled_restores_P2", "true");

    set_region_likelihood_rows_cache_diagnostic(2);
    clear_region_likelihood_rows_tls();
    let _ = diagnose(&pred_merged[0]);
    let pred_ck_inner = take_last_site_score_inner_trace().expect("pred ck");
    let content = diagnose(&cache_target[0]).expect("content key");
    let ck_inner = take_last_site_score_inner_trace().expect("ck inner");
    let pred_ck_lookup = pred_ck_inner.lookup.as_ref().expect("pred ck lookup");
    let ck_lookup = ck_inner.lookup.as_ref().expect("ck lookup");
    assert_ne!(
        (pred_ck_lookup.input_ptr, pred_ck_lookup.input_len),
        (ck_lookup.input_ptr, ck_lookup.input_len)
    );
    assert!(!ck_lookup.invalid_alias);
    assert_eq!(ck_inner.dense_hash, 0x7e3dbd79bc98a112);
    assert_eq!(content.genotype.format.ad_as_i32(), vec![42, 5]);
    assert_eq!(content.genotype.format.pl_as_i32(), vec![84, 0, 1738]);
    kv("variant_B_content_key_restores_P2", "true");
    kv(
        "variant_C_independent_storage",
        "already owned Vec<ReadLikelihoodRow>; extra copy cannot prevent key alias",
    );
    set_region_likelihood_rows_cache_diagnostic(0);
    clear_region_likelihood_rows_tls();
    kv(
        "decision_case",
        "B — legitimate intra-event reuse exists; pointer identity is not semantic",
    );
    kv(
        "minimal_semantic_contract",
        "HIT valid iff cached dense rows were rebuilt from the same sparse-cell population; pointer identity of a dropped subset Vec is not that population",
    );

    let closed_218_outcome = HaplotypeCallerEngine::call_region(
        covering_218,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call 218")
    .expect("outcome");
    assert_eq!(closed_218_outcome.assembly.haplotypes.len(), 30);
    let closed_218 = closed_218_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_218)
                && c.event.ref_allele == "G"
                && c.event.alt_allele == "T"
        })
        .expect("6R.218 G/T");
    assert_eq!(closed_218.genotype.format.pl_as_i32(), vec![69, 0, 2140]);
    let emitted_218 = try_emit_call_region_variants(
        covering_218,
        &closed_218_outcome,
        "SAMPLE",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("emit 218");
    let sib = emitted_218
        .iter()
        .find(|r| r.position == SIB && r.reference == "G")
        .expect("sib G/A");
    assert_eq!(
        sib.samples[0].pl.as_deref(),
        Some([645, 0, 1147].as_slice())
    );
    assert!((sib.quality.expect("sib QUAL") - 637.64).abs() < 0.01);
    kv("closed_6r218_pl", "69,0,2140");
    kv("closed_6r218_hap_n", "30");

    let tg_regions = flatten_assembly_regions(&walk(TG_INTERVAL, &p12_bam));
    let dp_regions = flatten_assembly_regions(&walk(DP_INTERVAL, &p12_bam));
    let mid_regions = flatten_assembly_regions(&walk(MID_INTERVAL, &p12_bam));
    let mid_b_regions = flatten_assembly_regions(&walk(MID_B_INTERVAL, &p12_bam));
    let post_regions = flatten_assembly_regions(&walk(P12_POST_INTERVAL, &p12_bam));
    let het_regions = flatten_assembly_regions(&walk(HET_TAIL_INTERVAL, &p12_bam));
    let covering_tg = tg_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_TG
                && r.end.get() >= CLOSED_TG
        })
        .expect("tg");
    let covering_ca = tg_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_CA
                && r.end.get() >= CLOSED_CA
        })
        .expect("ca");
    let covering_hom = mid_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_HOM_ALT
                && r.end.get() >= CLOSED_HOM_ALT
        })
        .expect("hom");
    let covering_one = mid_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_ONE_READ
                && r.end.get() >= CLOSED_ONE_READ
        })
        .expect("one");
    let covering_mid_b = mid_b_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_MID_B
                && r.end.get() >= CLOSED_MID_B
        })
        .expect("mid-B");
    let covering_post = post_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_POST
                && r.end.get() >= CLOSED_POST
        })
        .expect("post");
    let covering_het = het_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_WEAK
                && r.end.get() >= CLOSED_WEAK
        })
        .expect("het-tail");
    let covering_dp = dp_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_DP
                && r.end.get() >= CLOSED_DP
        })
        .expect("dp");
    let covering_ac = dp_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_AC
                && r.end.get() >= CLOSED_AC
        })
        .expect("ac");
    let call_one = |region: &gatk_haplotypecaller::AssemblyRegion| {
        HaplotypeCallerEngine::call_region(
            region,
            &dict,
            &ref_fasta,
            &CallRegionArgs::strict_java(),
        )
        .expect("call")
        .expect("outcome")
    };
    let tg_outcome = call_one(covering_tg);
    let ca_outcome = if std::ptr::eq(covering_tg, covering_ca) {
        tg_outcome.clone()
    } else {
        call_one(covering_ca)
    };
    let hom_outcome = call_one(covering_hom);
    let one_outcome = if std::ptr::eq(covering_hom, covering_one) {
        hom_outcome.clone()
    } else {
        call_one(covering_one)
    };
    let mid_b_outcome = call_one(covering_mid_b);
    let post_outcome = call_one(covering_post);
    let het_outcome = call_one(covering_het);
    let dp_outcome = call_one(covering_dp);
    let ac_outcome = if std::ptr::eq(covering_dp, covering_ac) {
        dp_outcome.clone()
    } else {
        call_one(covering_ac)
    };
    let find = |o: &gatk_haplotypecaller::CallRegionOutcome, pos: u64, r: &str, a: &str| {
        o.genotyped_calls
            .iter()
            .find(|c| {
                c.event.start_1based == GenomePosition::new_1based(pos)
                    && c.event.ref_allele == r
                    && c.event.alt_allele == a
            })
            .unwrap_or_else(|| panic!("{pos} {r}/{a}"))
            .clone()
    };
    let tg = find(&tg_outcome, CLOSED_TG, "T", "G");
    let ttc = find(&tg_outcome, CLOSED_TTC, "TTC", "T");
    let ct = find(&tg_outcome, CLOSED_CT, "CT", "C");
    let ca = find(&ca_outcome, CLOSED_CA, "C", "A");
    let hom = find(&hom_outcome, CLOSED_HOM_ALT, "A", "T");
    let del = find(&hom_outcome, CLOSED_DEL, "G", "A");
    let one = find(&one_outcome, CLOSED_ONE_READ, "C", "A");
    let mid_b = find(&mid_b_outcome, CLOSED_MID_B, "C", "A");
    let post = find(&post_outcome, CLOSED_POST, "C", "T");
    let het = find(&het_outcome, CLOSED_HET, "C", "T");
    let sib_p12 = find(&het_outcome, CLOSED_SIB, "G", "A");
    let weak = find(&het_outcome, CLOSED_WEAK, "C", "T");
    let dp = find(&dp_outcome, CLOSED_DP, "G", "T");
    let ag = find(&dp_outcome, CLOSED_AG, "A", "G");
    let ac = find(&ac_outcome, CLOSED_AC, "A", "C");
    assert_eq!(unique_indices(&tg.annotation_likelihoods).len(), 1);
    assert_eq!(ttc.genotype.format.pl_as_i32(), vec![45, 3, 0]);
    assert_eq!(ag.genotype.format.pl_as_i32(), vec![90, 6, 0]);
    assert_eq!(ac.genotype.format.pl_as_i32(), vec![130, 9, 0]);
    assert_eq!(dp.genotype.format.dp.as_i32(), 2);
    let _ = (ct, del);
    assert_eq!(unique_indices(&tg.annotation_likelihoods).len(), 1);
    assert_eq!(unique_indices(&ca.annotation_likelihoods).len(), 6);
    assert_eq!(unique_indices(&hom.annotation_likelihoods).len(), 2);
    assert_eq!(unique_indices(&one.annotation_likelihoods).len(), 1);
    assert_eq!(unique_indices(&mid_b.annotation_likelihoods).len(), 2);
    assert_eq!(unique_indices(&post.annotation_likelihoods).len(), 1);
    assert_eq!(unique_indices(&het.annotation_likelihoods).len(), 3);
    assert_eq!(unique_indices(&sib_p12.annotation_likelihoods).len(), 3);
    assert_eq!(unique_indices(&weak.annotation_likelihoods).len(), 3);
    assert_eq!(
        coverage_evidence_count(
            &het_outcome.genotyping_reads,
            &weak.annotation_likelihoods,
            CLOSED_WEAK,
            CLOSED_WEAK,
            MARGIN,
        ),
        3
    );
    assert_eq!(weak.genotype.format.pl_as_i32(), vec![55, 0, 21]);
    kv(
        "closed_controls",
        "6R.218 PL 69,0,2140 hap_n=30; P12 6R.174–6R.214; sib 29455019",
    );
}
