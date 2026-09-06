//! 6R.107 live dump: L9 post-emit-fail overwrite at HOLDOUT_6R53 `20:29455388 C/T`.
//! Skipped unless `HOLDOUT_6R107=1`.
//!
//! ```text
//! HOLDOUT_6R107=1 cargo test -p gatk-haplotypecaller --test holdout_6r107_l9_overwrite -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::VariationEvent;
use gatk_haplotypecaller::hc_genotyping_engine::{
    l9_may_overwrite_pairhmm_gls_after_emit_fail, SparsePlShape, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use serde_json::json;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const COVERING: (u64, u64) = (29_455_300, 29_455_559);
const TARGET: u64 = 29_455_388;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

#[test]
fn holdout_6r107_l9_overwrite_dump() {
    if std::env::var("HOLDOUT_6R107").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R107=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    assert!(ref_fasta.is_file(), "missing {}", ref_fasta.display());
    assert!(bam.is_file(), "missing {}", bam.display());

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
            ) && r.start.get() == COVERING.0
                && r.end.get() == COVERING.1
        })
        .expect("covering ActiveFull");
    let args = CallRegionArgs::strict_java();
    let outcome = HaplotypeCallerEngine::call_region(covering, &dict, &ref_fasta, &args)
        .expect("call")
        .expect("ActiveFull outcome");

    let event = VariationEvent::from_alleles("20", TARGET, "C", "T");
    let live = outcome.genotyped_calls.iter().find(|c| {
        c.event.start_1based.get() == TARGET
            && c.event.ref_allele == "C"
            && c.event.alt_allele == "T"
    });
    let pileup_ref = 44i32;
    let pileup_alt = 4i32;
    let overwrite = l9_may_overwrite_pairhmm_gls_after_emit_fail(&event, pileup_ref, pileup_alt);
    let emitted =
        try_emit_call_region_variants(covering, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let rust_emitted = emitted.iter().any(|r| r.position == TARGET);

    let doc = json!({
        "holdout": "20:29455388 C/T",
        "original_calculator_pl": [0, 6, 1780],
        "pileup_ref": pileup_ref,
        "pileup_alt": pileup_alt,
        "sparse_from_pileup": format!("{:?}", SparsePlShape::from_pileup_depths(pileup_ref, pileup_alt)),
        "l9_may_overwrite_pairhmm_after_emit_fail": overwrite,
        "live_genotyped": live.is_some(),
        "live_pl": live.map(|c| c.genotype.format.pl_as_i32()),
        "rust_emitted": rust_emitted,
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());

    assert!(!overwrite);
    assert!(
        live.is_none(),
        "calculator hom-ref must not be replaced with SparsePlShape"
    );
    assert!(!rust_emitted);
}
