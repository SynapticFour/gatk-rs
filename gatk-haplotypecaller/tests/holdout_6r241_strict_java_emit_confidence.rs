//! 6R.241 live: production strict-Java omit of covering `20:29455314 G>C`
//! at Java `stand-call-conf=30`. Same GLs still emit at threshold 10.
//! K=128 unchanged. Skipped unless `HOLDOUT_6R241=1`.
//!
//! ```text
//! HOLDOUT_6R241=1 cargo test -p gatk-haplotypecaller --test holdout_6r241_strict_java_emit_confidence -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::emit_gates::java_emit_would_pass;
use gatk_haplotypecaller::genotyping::best_pl_index;
use gatk_haplotypecaller::hc_genotyping_engine::{
    filter_genotyped_calls_for_strict_java_emit, java_emit_af_decision, HcGenotypingConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
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
const MATCHED: u64 = 29_455_379;
const COVERING: (u64, u64) = (29_455_300, 29_455_559);

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R241\t{key}\t{}", value.as_ref());
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

fn is_gc(c: &gatk_haplotypecaller::hc_genotyping_engine::GenotypedSiteCall) -> bool {
    c.event.start_1based.get() == TARGET && c.event.ref_allele == "G" && c.event.alt_allele == "C"
}

fn hit_gc(recs: &[gatk_core::io::vcf::VcfRecord]) -> bool {
    recs.iter()
        .any(|r| r.position == TARGET && r.reference == "G" && r.alternate.iter().any(|a| a == "C"))
}

#[test]
fn holdout_6r241_strict_java_emit_confidence() {
    if std::env::var("HOLDOUT_6R241").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R241=1");
        return;
    }
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128, "do not raise K");
    assert_eq!(
        HcGenotypingConfig::strict_java().stand_emit_confidence,
        DEFAULT_STAND_CALL_CONF
    );
    let root = repo_root();
    let java_vcf = root.join(JAVA_VCF_REL);
    assert!(java_vcf.is_file(), "missing {}", java_vcf.display());
    assert!(!java_vcf_has_gc(&java_vcf), "Java covering VCF omits G>C");

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

    let prod = HaplotypeCallerEngine::call_region(
        covering,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("prod")
    .expect("Some");
    assert!(prod
        .assembly
        .variation_events()
        .iter()
        .any(|e| { e.start_1based.get() == TARGET && e.ref_allele == "G" && e.alt_allele == "C" }));
    assert!(!prod.genotyped_calls.iter().any(is_gc));
    let prod_recs =
        try_emit_call_region_variants(covering, &prod, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    assert!(!hit_gc(&prod_recs));
    let ga = prod_recs
        .iter()
        .find(|r| {
            r.position == MATCHED && r.reference == "G" && r.alternate.iter().any(|a| a == "A")
        })
        .expect("G/A remains");
    let s = ga.samples.first().expect("sample");
    assert_eq!(s.ad.as_deref(), Some(&[42u32, 5][..]));
    assert_eq!(s.pl.as_deref(), Some(&[84u32, 0, 1738][..]));

    let mut args10 = CallRegionArgs::strict_java();
    args10.genotyping.stand_emit_confidence = DEFAULT_STAND_EMIT_CONFIDENCE;
    let at10 = HaplotypeCallerEngine::call_region(covering, &dict, &ref_fasta, &args10)
        .expect("at10")
        .expect("Some");
    let call = at10.genotyped_calls.iter().find(|c| is_gc(c)).expect("GLs");
    assert_eq!(call.genotype.format.pl_as_i32(), vec![21, 0, 1461]);
    assert_eq!(best_pl_index(&call.genotype.format.pl), 1);
    let gls = &call.genotype.genotype_log10_likelihoods;
    let af10 = java_emit_af_decision(gls, DEFAULT_STAND_EMIT_CONFIDENCE).expect("af10");
    let af30 = java_emit_af_decision(gls, DEFAULT_STAND_CALL_CONF).expect("af30");
    assert!(af10.passes_emit);
    assert!(!af30.passes_emit);
    assert!((af10.phred_scaled - 13.63).abs() < 0.02);
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
        try_emit_call_region_variants(covering, &at10, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    assert!(hit_gc(&recs10));

    let mut drop = vec![call.clone()];
    filter_genotyped_calls_for_strict_java_emit(
        &mut drop,
        &covering.reads,
        &at10.assembly,
        &HcGenotypingConfig::strict_java(),
    )
    .expect("filter");
    assert!(!drop.iter().any(is_gc));
    kv("classification", "EMISSION_PREDICATE_DIVERGENCE");
}
