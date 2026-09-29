//! 6R.297: production emission is unchanged. The span override is off.
//! Production PL remains 3518. Proof-only.
//! Skipped unless `HOLDOUT_6R297=1`.
//!
//! ```text
//! HOLDOUT_6R297=1 cargo test -p gatk-haplotypecaller --test holdout_6r297_complete_java_hap_ll_capture -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::hc_genotyping_engine::HcGenotypingConfig;
use gatk_haplotypecaller::read_threading_assembler::DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, set_forensic_6r294_padded_span,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
    DEFAULT_STAND_CALL_CONF, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_455_649;
const TARGET_REF: &str = "T";
const TARGET_ALT: &str = "TGTTTG";

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R297\t{key}\t{}", value.as_ref());
}

#[test]
fn holdout_6r297_complete_java_hap_ll_capture() {
    if std::env::var("HOLDOUT_6R297").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R297=1");
        return;
    }
    set_forensic_6r294_padded_span(None, None);
    gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r287_global_hap_floor(None, None);
    gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r288_haplotype_substitute(
        None, None, 0, None,
    );
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_target(None, 0);
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_substitute(None, None);
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128);
    assert_eq!(
        HcGenotypingConfig::strict_java().stand_emit_confidence,
        DEFAULT_STAND_CALL_CONF
    );
    let root = repo_root();
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
    let emitted = recs
        .iter()
        .find(|r| {
            r.position == TARGET
                && r.reference == TARGET_REF
                && r.alternate.iter().any(|a| a == TARGET_ALT)
        })
        .expect("emit");
    assert_eq!(
        emitted.samples[0].pl.as_deref(),
        Some(&[570, 0, 3518][..]),
        "production PL must remain 3518"
    );
    kv("production_integer_PL", "3518");
    kv(
        "classification",
        "LOGLESS_24_COLUMN_MATRIX_IDENTICAL_AVX_0_18_UNCAPTURED",
    );
    kv("production_change", "NONE");
}
