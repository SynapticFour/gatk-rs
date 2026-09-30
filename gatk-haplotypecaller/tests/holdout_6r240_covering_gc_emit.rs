//! 6R.240 live: covering `20:29455314 G>C` is in Java and Rust EventMaps
//! after 6R.239; Java covering VCF omits it. This holdout forces
//! `stand_emit=10` so the 10-vs-30 emit predicate stays observable;
//! production strict-Java 30 is 6R.241. K=128 unchanged.
//! Skipped unless `HOLDOUT_6R240=1`.
//!
//! ```text
//! HOLDOUT_6R240=1 cargo test -p gatk-haplotypecaller --test holdout_6r240_covering_gc_emit -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::emit_gates::java_emit_would_pass;
use gatk_haplotypecaller::hc_genotyping_engine::{
    java_emit_af_decision, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use gatk_haplotypecaller::read_threading_assembler::DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_CALL_CONF,
};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/java.vcf";
const TARGET: u64 = 29_455_314;
const COVERING: (u64, u64) = (29_455_300, 29_455_559);

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R240\t{key}\t{}", value.as_ref());
}

fn java_vcf_has_gc(path: &Path) -> bool {
    for line in std::fs::read_to_string(path).unwrap_or_default().lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<_> = line.split('\t').collect();
        if f.len() < 5 || f[0] != "20" {
            continue;
        }
        if f[1].parse::<u64>().ok() == Some(TARGET) && f[3] == "G" && f[4] == "C" {
            return true;
        }
    }
    false
}

#[test]
fn holdout_6r240_covering_gc_emit() {
    if std::env::var("HOLDOUT_6R240").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R240=1");
        return;
    }
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128, "do not raise K");
    let root = repo_root();
    let java_vcf = root.join(JAVA_VCF_REL);
    assert!(java_vcf.is_file(), "missing {}", java_vcf.display());
    assert!(
        !java_vcf_has_gc(&java_vcf),
        "pinned Java covering VCF omits 20:29455314 G>C"
    );
    kv("java_covering_vcf", "OMIT 20:29455314 G>C");

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
    let covering = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() == COVERING.0
                && r.end.get() == COVERING.1
        })
        .expect("covering");
    let mut args10 = CallRegionArgs::strict_java();
    args10.genotyping.stand_emit_confidence = DEFAULT_STAND_EMIT_CONFIDENCE;
    let outcome = HaplotypeCallerEngine::call_region(covering, &dict, &ref_fasta, &args10)
        .expect("call")
        .expect("Some");
    let has_event = outcome
        .assembly
        .variation_events()
        .iter()
        .any(|e| e.start_1based.get() == TARGET && e.ref_allele == "G" && e.alt_allele == "C");
    kv("rust_eventmap_gc", has_event.to_string());
    assert!(has_event, "Rust EventMap has G>C after 6R.239");

    let call = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based.get() == TARGET
                && c.event.ref_allele == "G"
                && c.event.alt_allele == "C"
        })
        .expect("Rust genotypes G>C");
    let gls = &call.genotype.genotype_log10_likelihoods;
    let af10 = java_emit_af_decision(gls, DEFAULT_STAND_EMIT_CONFIDENCE).expect("af10");
    let af30 = java_emit_af_decision(gls, DEFAULT_STAND_CALL_CONF).expect("af30");
    kv(
        "af",
        format!(
            "stand10_pass={} stand30_pass={} stand10_phred={:.2} stand30_mono={}",
            af10.passes_emit, af30.passes_emit, af10.phred_scaled, af30.site_is_monomorphic
        ),
    );
    assert!(af10.passes_emit);
    assert!(!af30.passes_emit);
    assert!(java_emit_would_pass(
        &call.event,
        gls,
        &call.genotype.format,
        DEFAULT_STAND_EMIT_CONFIDENCE,
        &[],
    )
    .unwrap());
    assert!(!java_emit_would_pass(
        &call.event,
        gls,
        &call.genotype.format,
        DEFAULT_STAND_CALL_CONF,
        &[],
    )
    .unwrap());

    let recs10 =
        try_emit_call_region_variants(covering, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let recs30 =
        try_emit_call_region_variants(covering, &outcome, "SAMPLE", DEFAULT_STAND_CALL_CONF)
            .unwrap_or_default();
    let hit = |recs: &[gatk_core::io::vcf::VcfRecord]| {
        recs.iter().any(|r| {
            r.position == TARGET && r.reference == "G" && r.alternate.iter().any(|a| a == "C")
        })
    };
    kv(
        "emit",
        format!("stand10={} stand30={}", hit(&recs10), hit(&recs30)),
    );
    assert!(hit(&recs10));
    assert!(!hit(&recs30));
    kv("classification", "EMISSION_PREDICATE_DIVERGENCE");
    kv("production_change", "NONE");
}
