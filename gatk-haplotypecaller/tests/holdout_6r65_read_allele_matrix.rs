//! 6R.65 forensic: genotyping read set / QNAME collapse at `20:29456344`.
//!
//! Skipped unless `HOLDOUT_6R65=1`.
//!
//! ```text
//! HOLDOUT_6R65=1 cargo test -p gatk-haplotypecaller --test holdout_6r65_read_allele_matrix -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, take_colocated_merge_numerics,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use serde_json::json;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const POS_SNP: u64 = 29_456_344;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

#[test]
fn holdout_6r65_read_allele_matrix_29456344() {
    if std::env::var("HOLDOUT_6R65").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R65=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    assert!(ref_fasta.is_file() && bam.is_file());

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
    let region = covering[0];
    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");

    let emitted =
        try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let near: Vec<_> = emitted
        .iter()
        .filter(|r| r.position >= POS_SNP.saturating_sub(5) && r.position <= POS_SNP + 5)
        .map(|r| {
            json!({
                "pos": r.position,
                "ref": r.reference,
                "alt": r.alternate,
                "gt": r.samples.get(0).and_then(|s| s.gt.as_ref().map(|g| g.alleles.clone())),
                "ad": r.samples.get(0).and_then(|s| s.ad.clone()),
                "pl": r.samples.get(0).and_then(|s| s.pl.clone()),
            })
        })
        .collect();
    let vcf_tc = emitted.iter().find(|r| {
        r.position == POS_SNP && r.reference == "T" && r.alternate.iter().any(|a| a == "C")
    });

    let live = take_colocated_merge_numerics();
    let numerics = live.iter().find(|n| n.loc == POS_SNP).cloned();
    let merge_fired = numerics.is_some();

    let doc = json!({
        "locus": "20:29456344 T/C",
        "merge_fired": merge_fired,
        "vcf_t_c": vcf_tc.map(|vcf| json!({
            "ref": vcf.reference,
            "alt": vcf.alternate,
            "gt": vcf.samples[0].gt.as_ref().map(|g| g.alleles.clone()),
            "ad": vcf.samples[0].ad.clone(),
            "pl": vcf.samples[0].pl.clone(),
        })),
        "vcf_near_pm5": near,
        "java_oracle_emitted_pl": [542, 0, 1353],
        "ladder_a": numerics.as_ref().map(|n| json!({
            "n_pairhmm_reads": n.n_pairhmm_reads,
            "n_overlap_before_qname_dedupe": n.n_overlap_before_qname_dedupe,
            "n_overlap_unique_qnames": n.n_overlap_unique_qnames,
            "n_qnames_with_multiple_overlapping_reads": n.n_qnames_with_multiple_overlapping_reads,
            "n_reads_after_qname_dedupe": n.n_reads,
            "n_haps": n.n_haps,
            "n_haps_with_multiple_events_at_loc": n.n_haps_with_multiple_events_at_loc,
            "pool_size_sum": n.pool_sizes.iter().sum::<usize>(),
        })),
        "ladder_b": numerics.as_ref().map(|n| json!({
            "rust_exclusive_pool_sizes": n.pool_sizes,
            "java_style_pool_sizes": n.java_style_pool_sizes,
            "hap_event_signatures_at_loc": n.hap_event_signatures_at_loc,
            "n_cache_events": n.n_cache_events,
            "nearest_event_start_below": n.nearest_event_start_below,
            "nearest_event_start_above": n.nearest_event_start_above,
        })),
        "ladder_c": numerics.as_ref().map(|n| json!({
            "n_allele_floor_clips": n.n_allele_floor_clips,
            "merged_pl": n.merged_pl,
            "subset_pl": n.subset_pl,
        })),
        "note": "6R.66: empty-EventMap pad-slice removed. Merge NotApplicable (no live numerics) means alt pools stayed empty — Java-compatible EventMap mapping.",
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());

    if let Some(n) = numerics {
        assert_ne!(
            n.pool_sizes,
            vec![40, 6, 21],
            "6R.66: pad-slice pools 40/6/21 must not remain"
        );
        if !n.java_style_pool_sizes.is_empty() {
            assert_eq!(
                n.pool_sizes, n.java_style_pool_sizes,
                "production exclusive pools must match Java EventMap-only mapper"
            );
        }
    }
}
