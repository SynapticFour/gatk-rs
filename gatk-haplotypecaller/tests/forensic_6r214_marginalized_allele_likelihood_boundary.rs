//! 6R.214: first causal arrow in the marginalized allele-likelihood
//! object at `20:29455015 G/T`.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE. Pre-marginalization haplotype population already
//! differs (Java 30 vs Rust 84). Do not retune QUAL/GQ/QD or compensate in
//! `marginalize`.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r214_marginalized_allele_likelihood_boundary -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{
    build_per_haplotype_variation_events, variation_events_at_position_from_cache,
};
use gatk_haplotypecaller::haplotype::Haplotype;
use gatk_haplotypecaller::hc_allele_mapping::create_allele_mapper_with_events;
use gatk_haplotypecaller::hc_genotyping_engine::{
    alt_hap_indices_for_genotype_marginalization, biallelic_genotype_log10_likelihoods_gatk,
    java_alignment_read_overlaps_interval, marginalize_rows_to_biallelic_alleles,
    ref_hap_indices_for_genotype_marginalization, region_likelihoods_to_rows, HcGenotypingConfig,
    DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, ReadLikelihoodRow, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_455_015;
const SIB: u64 = 29_455_019;
const JAVA_TRIM_START: u64 = 29_454_995;
const JAVA_TRIM_END: u64 = 29_455_124;
const TARGET_REF: &str = "G";
const TARGET_ALT: &str = "T";

/// Java `genotype-emit-at-loc` / `ll-input-at-loc` FNV-1a of haplotype bases.
/// Assembly output is 30 trimmed 130 bp haplotypes (`alignStart=595`) *before*
/// PairHMM, so GKL vs LOGLESS does not change this population.
const JAVA_HAP_HASHES: &[&str] = &[
    "a9a6a6a3c7e334fe",
    "17b253b3caefbc24",
    "8cdf90fd00055da9",
    "94152d17f55d2c0f",
    "504af452f2b4a369",
    "421f004a656474e4",
    "f61271afcbc2fb1d",
    "dd22d427216a7c13",
    "a075ca46f85f60ba",
    "66751242bb88dd8c",
    "f5d1e8d60cfb9666",
    "1ec62235b5a6f477",
    "ba566a88f4a37b49",
    "f9fb0ec02bc3a8ef",
    "e24a3d9b7270ab1e",
    "a288c9e831abc8c3",
    "65773746fc8fc8cd",
    "80b1e7bc60742470",
    "8a328239b600c7bd",
    "cf7d9939fecc0fb3",
    "fad3d8d2e57bba5a",
    "a1c7814baa41b868",
    "04cd4d695bcfcef2",
    "a911d3c74c6ac1ab",
    "7722ec426dd6775e",
    "8183d80853176284",
    "405bde68757e1909",
    "7c0149c997f1c2d1",
    "85fb4ea1d9d0a2b7",
    "fdfee8cd394ccda6",
];

const JAVA_REF_HASHES: &[&str] = &[
    "a9a6a6a3c7e334fe",
    "17b253b3caefbc24",
    "8cdf90fd00055da9",
    "94152d17f55d2c0f",
    "504af452f2b4a369",
    "421f004a656474e4",
    "f61271afcbc2fb1d",
    "dd22d427216a7c13",
    "a075ca46f85f60ba",
    "66751242bb88dd8c",
    "f5d1e8d60cfb9666",
    "1ec62235b5a6f477",
    "ba566a88f4a37b49",
    "f9fb0ec02bc3a8ef",
    "e24a3d9b7270ab1e",
    "a288c9e831abc8c3",
    "65773746fc8fc8cd",
    "80b1e7bc60742470",
];

const JAVA_ALT_HASHES: &[&str] = &[
    "8a328239b600c7bd",
    "cf7d9939fecc0fb3",
    "fad3d8d2e57bba5a",
    "a1c7814baa41b868",
    "04cd4d695bcfcef2",
    "a911d3c74c6ac1ab",
    "7722ec426dd6775e",
    "8183d80853176284",
    "405bde68757e1909",
    "7c0149c997f1c2d1",
    "85fb4ea1d9d0a2b7",
    "fdfee8cd394ccda6",
];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R214\t{key}\t{}", value.as_ref());
}

fn fnv1a64_hex(bases: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bases {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{h:016x}")
}

fn window_bases(h: &Haplotype, win_start: u64, win_end: u64) -> Option<Vec<u8>> {
    let gl = h.genome_loc?;
    let gs = gl.start_1based();
    if win_start < gs || win_end < win_start {
        return None;
    }
    let off = (win_start - gs) as usize;
    let len = (win_end - win_start + 1) as usize;
    h.bases.get(off..off.checked_add(len)?).map(|s| s.to_vec())
}

fn allele_of(idx: usize, ref_set: &BTreeSet<usize>, alt_set: &BTreeSet<usize>) -> &'static str {
    if alt_set.contains(&idx) {
        "T"
    } else if ref_set.contains(&idx) {
        "G*"
    } else {
        "unmapped"
    }
}

#[test]
fn forensic_6r214_java_pre_marginalize_object_pin() {
    kv("java_pin", JAVA_PIN);
    kv("java_pairhmm_input_hap_n", "30");
    kv("java_assign_gl_hap_n", "30");
    kv("java_hap_len", "130");
    kv("java_hap_align_start", "595");
    kv("java_trim_span", "20:29454995-29455124");
    kv("java_hap_ll_evidence", "88");
    kv("java_retain_evidence", "65");
    kv("java_mapper_keys", "G*,T");
    kv("java_mapper_ref_n", "18");
    kv("java_mapper_alt_n", "12");
    kv("java_pl", "69,0,2140");
    kv("java_log10_gl", "-6.9,0,-214");
    kv("java_qual", "61.63958660169907");
    kv("java_ref_hash", "e24a3d9b7270ab1e");
    assert_eq!(JAVA_HAP_HASHES.len(), 30);
    assert_eq!(JAVA_REF_HASHES.len(), 18);
    assert_eq!(JAVA_ALT_HASHES.len(), 12);
    assert_eq!(
        JAVA_REF_HASHES
            .iter()
            .chain(JAVA_ALT_HASHES.iter())
            .copied()
            .collect::<BTreeSet<_>>(),
        JAVA_HAP_HASHES.iter().copied().collect()
    );
}

#[test]
fn forensic_6r214_live_haplotype_likelihood_input_boundary() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }
    kv("java_pin", JAVA_PIN);
    kv(
        "production_change",
        "NONE in 6R.214; closed upstream by 6R.218 GRAPH_STATE_DIVERGENCE",
    );
    kv(
        "first_divergent_arrow",
        "closed: haplotype population now matches Java 30 after 6R.218 cycle abort",
    );
    kv("classification", "CLOSED_AFTER_6R218");

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
    assert_eq!(covering.len(), 1, "one ActiveFull covers the target");
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
    let config = HcGenotypingConfig::strict_java();
    let emit_spanning = !config.disable_spanning_event_genotyping;
    let (full_ref, full_pad) = outcome.assembly.event_map_reference();
    let hap_events = build_per_haplotype_variation_events(
        haps,
        full_ref,
        full_pad,
        outcome.assembly.max_mnp_distance(),
        &region.contig,
    );
    let cache_only = variation_events_at_position_from_cache(&hap_events, TARGET, emit_spanning);
    assert_eq!(cache_only.len(), 1);
    assert_eq!(cache_only[0].ref_allele, TARGET_REF);
    assert_eq!(cache_only[0].alt_allele, TARGET_ALT);

    let apply_bases = outcome.assembly.apply_bases_shared();
    let apply_pad = haps
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.map(|g| g.start_1based()))
        .unwrap_or(full_pad);
    let mapping = create_allele_mapper_with_events(
        &cache_only[0],
        TARGET,
        haps,
        apply_pad,
        apply_bases.as_ref(),
        outcome.assembly.max_mnp_distance(),
        emit_spanning,
        Some(&hap_events),
    );
    let ref_set: BTreeSet<usize> = mapping
        .ref_haplotype_indices
        .iter()
        .map(|i| i.get())
        .collect();
    let alt_set: BTreeSet<usize> = mapping
        .alt_haplotype_indices
        .iter()
        .map(|i| i.get())
        .collect();
    let ref_hap = haps.iter().find(|h| h.is_reference).expect("ref hap");
    let prod_ref =
        ref_hap_indices_for_genotype_marginalization(&mapping, haps, &config, Some(&cache_only[0]));
    let prod_alt = alt_hap_indices_for_genotype_marginalization(
        &mapping,
        haps,
        &cache_only[0],
        ref_hap,
        apply_pad,
        apply_bases.as_ref(),
        outcome.assembly.max_mnp_distance(),
        &region.contig,
        &config,
    );

    let mut rust_full = BTreeSet::new();
    let mut rust_window = BTreeSet::new();
    let mut rust_len: BTreeMap<usize, usize> = BTreeMap::new();
    let mut rust_span: BTreeMap<String, usize> = BTreeMap::new();
    let mut window_map: BTreeMap<String, String> = BTreeMap::new();
    for (i, h) in haps.iter().enumerate() {
        let full = fnv1a64_hex(&h.bases);
        rust_full.insert(full.clone());
        *rust_len.entry(h.bases.len()).or_insert(0) += 1;
        let span = h
            .genome_loc
            .map(|g| format!("{}-{}", g.start_1based(), g.end_1based()))
            .unwrap_or_else(|| ".".to_string());
        *rust_span.entry(span).or_insert(0) += 1;
        let win_hash = window_bases(h, JAVA_TRIM_START, JAVA_TRIM_END)
            .map(|b| fnv1a64_hex(&b))
            .unwrap_or_else(|| format!("no_window:{full}"));
        rust_window.insert(win_hash.clone());
        window_map.insert(
            full,
            format!("{win_hash}\t{}", allele_of(i, &ref_set, &alt_set)),
        );
        kv(
            "rust_hap",
            format!(
                "{i}\t{}\tlen={}\tisRef={}\tmapped={}\twindow={}",
                fnv1a64_hex(&h.bases),
                h.bases.len(),
                h.is_reference,
                allele_of(i, &ref_set, &alt_set),
                win_hash
            ),
        );
    }

    let java: BTreeSet<&str> = JAVA_HAP_HASHES.iter().copied().collect();
    let rust_full_ref: BTreeSet<&str> = rust_full.iter().map(String::as_str).collect();
    let rust_win_ref: BTreeSet<&str> = rust_window.iter().map(String::as_str).collect();
    let common_full: Vec<_> = java.intersection(&rust_full_ref).copied().collect();
    let java_only_full: Vec<_> = java.difference(&rust_full_ref).copied().collect();
    let rust_only_full: Vec<_> = rust_full_ref.difference(&java).copied().collect();
    let common_win: Vec<_> = java.intersection(&rust_win_ref).copied().collect();
    let java_only_win: Vec<_> = java.difference(&rust_win_ref).copied().collect();
    let rust_only_win: Vec<_> = rust_win_ref.difference(&java).copied().collect();

    kv("java_hap_n", JAVA_HAP_HASHES.len().to_string());
    kv("rust_hap_n", haps.len().to_string());
    kv("rust_unique_full_hash_n", rust_full.len().to_string());
    kv("rust_unique_window_hash_n", rust_window.len().to_string());
    kv("common_full_hash_n", common_full.len().to_string());
    kv("java_only_full_hash_n", java_only_full.len().to_string());
    kv("rust_only_full_hash_n", rust_only_full.len().to_string());
    kv("common_window_hash_n", common_win.len().to_string());
    kv("java_only_window_hash_n", java_only_win.len().to_string());
    kv("rust_only_window_hash_n", rust_only_win.len().to_string());
    kv(
        "rust_hap_lens",
        rust_len
            .iter()
            .map(|(k, v)| format!("{k}={v}"))
            .collect::<Vec<_>>()
            .join(","),
    );
    kv(
        "rust_genome_spans",
        rust_span
            .iter()
            .map(|(k, v)| format!("{k}={v}"))
            .collect::<Vec<_>>()
            .join(","),
    );
    kv(
        "mapper_ref_n",
        mapping.ref_haplotype_indices.len().to_string(),
    );
    kv(
        "mapper_alt_n",
        mapping.alt_haplotype_indices.len().to_string(),
    );
    kv("prod_ref_pool_n", prod_ref.len().to_string());
    kv("prod_alt_pool_n", prod_alt.len().to_string());
    kv("java_mapper_ref_n", "18");
    kv("java_mapper_alt_n", "12");

    let rows = region_likelihoods_to_rows(&outcome.read_likelihoods, haps.len());
    let overlap_idx: BTreeSet<usize> = outcome
        .genotyping_reads
        .iter()
        .enumerate()
        .filter(|(_, rec)| {
            java_alignment_read_overlaps_interval(
                rec,
                TARGET,
                TARGET,
                DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN,
            )
        })
        .map(|(i, _)| i)
        .collect();
    let overlap_rows: Vec<ReadLikelihoodRow> = rows
        .iter()
        .filter(|r| overlap_idx.contains(&r.read_index))
        .cloned()
        .collect();
    let marg = marginalize_rows_to_biallelic_alleles(&overlap_rows, &prod_ref, &prod_alt);
    let gls = biallelic_genotype_log10_likelihoods_gatk(&marg, 0, 1);

    kv(
        "rust_genotyping_read_n",
        outcome.genotyping_reads.len().to_string(),
    );
    kv("rust_likelihood_row_n", rows.len().to_string());
    kv("rust_overlap_read_n", overlap_idx.len().to_string());
    kv("java_hap_ll_evidence", "88");
    kv("java_retain_evidence", "65");
    kv(
        "overlap_gls",
        gls.iter()
            .map(|x| format!("{x:.8}"))
            .collect::<Vec<_>>()
            .join(","),
    );

    let site = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based.get() == TARGET
                && c.event.ref_allele == TARGET_REF
                && c.event.alt_allele == TARGET_ALT
        })
        .expect("G/T");
    kv(
        "rust_pl",
        site.genotype
            .format
            .pl_as_i32()
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(","),
    );
    kv(
        "downstream_gls",
        site.genotype
            .genotype_log10_likelihoods
            .iter()
            .map(|x| format!("{x:.8}"))
            .collect::<Vec<_>>()
            .join(","),
    );

    assert_eq!(
        haps.len(),
        30,
        "6R.218 restored Java PairHMM/call_region hap_n"
    );
    assert_eq!(JAVA_HAP_HASHES.len(), 30);
    assert_eq!(
        common_win.len(),
        30,
        "every Java hap is present as a Rust window"
    );
    assert!(java_only_win.is_empty());
    assert!(
        rust_only_win.is_empty(),
        "Rust-only trim-window sequences must vanish with cyclic k=10 abort"
    );
    assert_eq!(mapping.ref_haplotype_indices.len(), 18);
    assert_eq!(mapping.alt_haplotype_indices.len(), 12);
    assert_eq!(prod_ref.len(), 18);
    assert_eq!(prod_alt.len(), 12);
    assert_eq!(site.genotype.format.pl_as_i32(), vec![69, 0, 2140]);

    let emitted =
        try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let sib = emitted
        .iter()
        .find(|r| r.position == SIB && r.reference == "G" && r.alternate.iter().any(|a| a == "A"))
        .expect("sib");
    assert_eq!(
        sib.samples[0].pl.as_deref(),
        Some([645, 0, 1147].as_slice())
    );
    assert!((sib.quality.expect("sib QUAL") - 637.64).abs() < 0.01);

    kv(
        "causal_proof",
        "6R.218: PairHMM/call_region hap_n=30 (18 REF/12 ALT) matches Java; PL 69,0,2140.",
    );
    kv(
        "not_investigated",
        "marginalize max vs log-sum-exp, post-marginalize -4.5 floor, GL calculator (Case C/D blocked by Case A)",
    );
    let _ = (
        window_map,
        overlap_rows,
        java_only_full,
        rust_only_full,
        java_only_win,
        rust_only_win,
        common_win,
    );
}
