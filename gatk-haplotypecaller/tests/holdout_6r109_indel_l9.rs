//! 6R.109 live dump: indel L9 must not overwrite valid hom-ref calculator GLs
//! at HOLDOUT_6R53 remainder `20:29455644 A/AT`. Skipped unless `HOLDOUT_6R109=1`.
//!
//! Adjacent `20:29455649` is recorded only as collateral; it is not this arrow.
//!
//! ```text
//! HOLDOUT_6R109=1 cargo test -p gatk-haplotypecaller --test holdout_6r109_indel_l9 -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::VariationEvent;
use gatk_haplotypecaller::hc_genotyping_engine::{
    l9_may_overwrite_pairhmm_gls_after_emit_fail, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use serde_json::json;
use std::fs;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/java.vcf";
const REGION: (u64, u64) = (29_455_560, 29_455_744);
const TARGET: u64 = 29_455_644;
const ADJACENT: u64 = 29_455_649;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn vcf_alleles_at(path: &Path, pos: u64) -> Vec<String> {
    let mut out = Vec::new();
    if !path.is_file() {
        return out;
    }
    for line in fs::read_to_string(path).unwrap_or_default().lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<_> = line.split('\t').collect();
        if f.len() < 5 || f[0] != "20" {
            continue;
        }
        if f[1].parse::<u64>().ok() == Some(pos) {
            out.push(format!("{}/{}", f[3], f[4]));
        }
    }
    out
}

#[test]
fn holdout_6r109_indel_l9_live_absent() {
    if std::env::var("HOLDOUT_6R109").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R109=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    let java_vcf = root.join(JAVA_VCF_REL);
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
            ) && r.start.get() == REGION.0
                && r.end.get() == REGION.1
        })
        .expect("ActiveFull 20:29455560-29455744");
    let args = CallRegionArgs::strict_java();
    let outcome = HaplotypeCallerEngine::call_region(covering, &dict, &ref_fasta, &args)
        .expect("call")
        .expect("ActiveFull outcome");

    let event = VariationEvent::from_alleles("20", TARGET, "A", "AT");
    let live_644 = outcome.genotyped_calls.iter().find(|c| {
        c.event.start_1based.get() == TARGET
            && c.event.ref_allele == "A"
            && c.event.alt_allele == "AT"
    });
    let live_649: Vec<_> = outcome
        .genotyped_calls
        .iter()
        .filter(|c| c.event.start_1based.get() == ADJACENT)
        .map(|c| {
            json!({
                "ref": c.event.ref_allele,
                "alt": c.event.alt_allele,
                "pl": c.genotype.format.pl_as_i32(),
                "ad": c.genotype.format.ad_as_i32(),
            })
        })
        .collect();
    let overwrite = l9_may_overwrite_pairhmm_gls_after_emit_fail(&event, 8, 4, true);
    let emitted =
        try_emit_call_region_variants(covering, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let rust_644 = emitted.iter().any(|r| {
        r.position == TARGET && r.reference == "A" && r.alternate.iter().any(|a| a == "AT")
    });
    let rust_649: Vec<_> = emitted
        .iter()
        .filter(|r| r.position == ADJACENT)
        .map(|r| format!("{}/{}", r.reference, r.alternate.join(",")))
        .collect();

    let java_644 = vcf_alleles_at(&java_vcf, TARGET);
    let java_649 = vcf_alleles_at(&java_vcf, ADJACENT);

    let doc = json!({
        "holdout": "20:29455644 A/AT",
        "active_full": REGION,
        "java_frozen_644": java_644,
        "java_frozen_649": java_649,
        "l9_may_overwrite_valid_homref_indel": overwrite,
        "live_genotyped_644": live_644.is_some(),
        "live_pl_644": live_644.map(|c| c.genotype.format.pl_as_i32()),
        "rust_emitted_644": rust_644,
        "live_genotyped_649": live_649,
        "rust_emitted_649": rust_649,
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());

    assert!(
        java_644.is_empty(),
        "frozen Java VCF must not contain 20:29455644"
    );
    assert!(!overwrite);
    assert!(
        live_644.is_none(),
        "valid hom-ref calculator must not be replaced with SparsePlShape"
    );
    assert!(!rust_644, "Rust must not emit 20:29455644 A/AT");
}
