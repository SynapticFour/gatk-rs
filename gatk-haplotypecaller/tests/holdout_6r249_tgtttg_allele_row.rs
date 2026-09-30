//! 6R.249 live: `20:29455649 T/TGTTTG` allele-row is max over five
//! EventMap-mapped TGTTTG haplotypes. Proof-only. Skipped unless
//! `HOLDOUT_6R249=1`.
//!
//! ```text
//! HOLDOUT_6R249=1 cargo test -p gatk-haplotypecaller --test holdout_6r249_tgtttg_allele_row -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{
    build_per_haplotype_variation_events, overlapping_events, remap_alt_onto_longer_ref,
};
use gatk_haplotypecaller::hc_genotyping_engine::{
    take_colocated_merge_numerics, HcGenotypingConfig,
};
use gatk_haplotypecaller::read_threading_assembler::DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, region_likelihoods_to_rows,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, GenomePosition, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
    DEFAULT_STAND_CALL_CONF, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/java.vcf";
const TARGET: u64 = 29_455_649;
const TARGET_REF: &str = "T";
const TARGET_ALT: &str = "TGTTTG";
const UNUSED_ALT: &str = "TTTG";
const CLOSED_GC: u64 = 29_455_314;
const EMPTY_POOL: f64 = -50.0;
const FLOOR_CAP: f64 = -4.5;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R249\t{key}\t{}", value.as_ref());
}

fn vcf_has(path: &Path, pos: u64, r: &str, a: &str) -> bool {
    for line in std::fs::read_to_string(path).unwrap_or_default().lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<_> = line.split('\t').collect();
        if f.len() < 5 || f[0] != "20" {
            continue;
        }
        if f[1].parse::<u64>().ok() == Some(pos) && f[3] == r && f[4] == a {
            return true;
        }
    }
    false
}

fn java_mapper_pools(
    n_haps: usize,
    hap_events: &gatk_haplotypecaller::event_map::PerHaplotypeVariationEvents,
    loc: u64,
    long_ref: &str,
    alts: &[String],
) -> Vec<Vec<usize>> {
    let loc_pos = GenomePosition::new_1based(loc);
    let mut pools: Vec<Vec<usize>> = vec![Vec::new(); 1 + alts.len()];
    for i in 0..n_haps {
        let spanning = overlapping_events(hap_events.events_for(i), loc);
        if spanning.is_empty() {
            pools[0].push(i);
            continue;
        }
        let mut in_alt = vec![false; alts.len()];
        let mut in_ref = false;
        for ev in &spanning {
            if ev.start_1based == loc_pos {
                if ev.ref_allele.len() == long_ref.len() {
                    if let Some(ai) = alts.iter().position(|a| a == &ev.alt_allele) {
                        in_alt[ai] = true;
                    }
                } else if ev.ref_allele.len() < long_ref.len() {
                    if let Some(remapped) =
                        remap_alt_onto_longer_ref(&ev.ref_allele, &ev.alt_allele, long_ref)
                    {
                        if let Some(ai) = alts.iter().position(|a| a == &remapped) {
                            in_alt[ai] = true;
                        }
                    }
                }
            } else {
                if let Some(ai) = alts.iter().position(|a| a == "*") {
                    in_alt[ai] = true;
                }
                break;
            }
        }
        if !in_ref && !in_alt.iter().any(|&b| b) {
            continue;
        }
        if in_ref {
            pools[0].push(i);
        }
        for (ai, hit) in in_alt.iter().enumerate() {
            if *hit {
                pools[ai + 1].push(i);
            }
        }
    }
    pools
}

fn pool_max(indices: &[usize], lls: &[f64]) -> f64 {
    let ll = indices
        .iter()
        .filter_map(|&i| lls.get(i).copied())
        .fold(f64::NEG_INFINITY, f64::max);
    if ll.is_finite() {
        ll
    } else {
        EMPTY_POOL
    }
}

#[test]
fn holdout_6r249_tgtttg_allele_row() {
    if std::env::var("HOLDOUT_6R249").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R249=1");
        return;
    }
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128);
    assert_eq!(
        HcGenotypingConfig::strict_java().stand_emit_confidence,
        DEFAULT_STAND_CALL_CONF
    );
    let root = repo_root();
    let java_vcf = root.join(JAVA_VCF_REL);
    assert!(vcf_has(&java_vcf, TARGET, TARGET_REF, TARGET_ALT));
    assert!(!vcf_has(&java_vcf, CLOSED_GC, "G", "C"));

    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
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
    let recs =
        try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    assert!(
        !recs.iter().any(|r| r.position == CLOSED_GC
            && r.reference == "G"
            && r.alternate.iter().any(|a| a == "C")),
        "covering G>C must stay omitted from this emit"
    );
    let snap = take_colocated_merge_numerics()
        .into_iter()
        .find(|s| s.loc == TARGET)
        .expect("merge snapshot");

    kv("target", "20:29455649 T/TGTTTG");
    kv("production_change", "NONE");
    kv("pool_sizes", format!("{:?}", snap.pool_sizes));
    kv(
        "java_style_pool_sizes",
        format!("{:?}", snap.java_style_pool_sizes),
    );
    assert_eq!(snap.assigned_gt, vec![0, 2]);
    assert_eq!(snap.remaining_keep_indices, vec![0, 2]);
    assert_eq!(snap.alts, [UNUSED_ALT, TARGET_ALT]);
    assert_eq!(snap.n_reads, 123);
    assert_eq!(snap.pool_sizes, snap.java_style_pool_sizes);
    assert_eq!(snap.pool_sizes[2], 5);
    assert_eq!(snap.n_allele_floor_clips, 0);
    assert_eq!(snap.subset_pl, vec![570, 0, 3518]);

    let haps = &outcome.assembly.haplotypes;
    let hap_events = build_per_haplotype_variation_events(
        haps,
        outcome.assembly.reference_bases(),
        outcome.assembly.padded_reference_start_1based(),
        outcome.assembly.max_mnp_distance(),
        "20",
    );
    let alts = vec![UNUSED_ALT.to_string(), TARGET_ALT.to_string()];
    let pools = java_mapper_pools(haps.len(), &hap_events, TARGET, "T", &alts);
    assert_eq!(pools[2].len(), 5);
    let loc_pos = GenomePosition::new_1based(TARGET);
    for &hi in &pools[2] {
        let remapped: Vec<String> = overlapping_events(hap_events.events_for(hi), TARGET)
            .iter()
            .filter(|e| e.start_1based == loc_pos)
            .filter_map(|e| {
                if e.ref_allele.len() == 1 {
                    Some(e.alt_allele.clone())
                } else {
                    remap_alt_onto_longer_ref(&e.ref_allele, &e.alt_allele, "T")
                }
            })
            .collect();
        assert!(
            remapped.iter().any(|a| a == TARGET_ALT),
            "pool member {hi} must carry EventMap TGTTTG"
        );
        assert!(!haps[hi].is_reference);
    }

    let retain: HashSet<usize> = snap.ad_row_read_index.iter().copied().collect();
    let subset: Vec<_> = outcome
        .read_likelihoods
        .iter()
        .filter(|c| retain.contains(&c.read_index.get()))
        .cloned()
        .collect();
    let hap_rows = region_likelihoods_to_rows(&subset, haps.len());
    let mut by_read = HashMap::new();
    for row in &hap_rows {
        by_read.entry(row.read_index).or_insert(row);
    }
    let mut n_unique = 0usize;
    let mut n_ties = 0usize;
    for (ri, &read_idx) in snap.ad_row_read_index.iter().enumerate() {
        let row = by_read.get(&read_idx).expect("row");
        let raw = pool_max(&pools[2], &row.haplotype_log10_likelihoods);
        let n_at = pools[2]
            .iter()
            .filter(|&&hi| row.haplotype_log10_likelihoods.get(hi).copied() == Some(raw))
            .count();
        if n_at == 1 {
            n_unique += 1;
        } else if n_at > 1 {
            n_ties += 1;
        }
        let pooled: Vec<f64> = pools
            .iter()
            .map(|p| pool_max(p, &row.haplotype_log10_likelihoods))
            .collect();
        let best = pooled.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let floor = best + FLOOR_CAP;
        let after = if pooled[2].is_finite() && pooled[2] < floor {
            floor
        } else {
            pooled[2]
        };
        assert!(
            (after - snap.ad_row_lls[ri][2]).abs() <= 1e-15,
            "reconstructed TGTTTG column must match ad_row"
        );
    }
    kv("n_unique_max", n_unique.to_string());
    kv("n_ties", n_ties.to_string());
    assert_eq!(n_unique + n_ties, 123);
    kv("pooling", "max over five EventMap TGTTTG haplotypes");
    kv("classification", "ALLELE_LIKELIHOOD_INPUT_DIVERGENCE");
    kv("production_change", "NONE");
}
