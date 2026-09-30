//! 6R.241: wire Java `standardConfidenceForCalling=30` into strict-Java emit.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! Same covering GLs as 6R.240 (`GT 0/1`, `PL 21,0,1461`, QUAL ~13.63).
//! Threshold 10 emits; production strict-Java 30 does not.
//! 6R.239 sentinel arithmetic is unchanged. K=128 is unchanged.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r241_strict_java_emit_confidence -- --nocapture --test-threads=1
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

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/java.vcf";
const TARGET: u64 = 29_455_314;
const TARGET_REF: &str = "G";
const TARGET_ALT: &str = "C";
const MATCHED: u64 = 29_455_379;
const COVERING: (u64, u64) = (29_455_300, 29_455_559);

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R241\t{key}\t{}", value.as_ref());
}

fn vcf_has_allele(path: &Path, pos: u64, r: &str, a: &str) -> bool {
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

fn rec_is(r: &gatk_core::io::vcf::VcfRecord, pos: u64, reference: &str, alt: &str) -> bool {
    r.position == pos && r.reference == reference && r.alternate.iter().any(|a| a == alt)
}

fn call_is_gc(c: &gatk_haplotypecaller::hc_genotyping_engine::GenotypedSiteCall) -> bool {
    c.event.start_1based.get() == TARGET
        && c.event.ref_allele == TARGET_REF
        && c.event.alt_allele == TARGET_ALT
}

fn args_at_stand(stand: f64) -> CallRegionArgs {
    let mut args = CallRegionArgs::strict_java();
    args.genotyping.stand_emit_confidence = stand;
    args
}

#[test]
fn forensic_6r241_strict_java_uses_java_calling_confidence() {
    kv("java_pin", JAVA_PIN);
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128, "do not raise K");
    assert_eq!(DEFAULT_STAND_EMIT_CONFIDENCE, 10.0);
    assert_eq!(DEFAULT_STAND_CALL_CONF, 30.0);
    let strict = HcGenotypingConfig::strict_java();
    assert_eq!(
        strict.stand_emit_confidence, DEFAULT_STAND_CALL_CONF,
        "strict-Java path must supply Java standardConfidenceForCalling=30"
    );

    let root = repo_root();
    let java_vcf = root.join(JAVA_VCF_REL);
    if !java_vcf.is_file() {
        eprintln!("skip: missing pinned Java covering VCF");
        return;
    }
    assert!(
        !vcf_has_allele(&java_vcf, TARGET, TARGET_REF, TARGET_ALT),
        "Java covering VCF omits 20:29455314 G>C"
    );
    kv("java_covering_vcf", "OMIT 20:29455314 G>C");

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
    let region = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() == COVERING.0
                && r.end.get() == COVERING.1
        })
        .expect("covering ActiveFull");

    let outcome30 = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call30")
    .expect("Some");
    let has_event = outcome30.assembly.variation_events().iter().any(|e| {
        e.start_1based.get() == TARGET && e.ref_allele == TARGET_REF && e.alt_allele == TARGET_ALT
    });
    kv("rust_eventmap_gc", has_event.to_string());
    assert!(has_event, "EventMap still has G>C after 6R.239");
    assert!(
        !outcome30.genotyped_calls.iter().any(call_is_gc),
        "production strict-Java calculateGenotypes omits G>C"
    );

    let recs30 =
        try_emit_call_region_variants(region, &outcome30, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    assert!(
        !recs30
            .iter()
            .any(|r| rec_is(r, TARGET, TARGET_REF, TARGET_ALT)),
        "production VCF emit omits G>C (emission decision only)"
    );
    let ga30 = recs30
        .iter()
        .find(|r| rec_is(r, MATCHED, "G", "A"))
        .expect("20:29455379 G/A remains emitted");
    let s30 = ga30.samples.first().expect("sample");
    kv(
        "matched_ga_format",
        format!(
            "GT={:?} AD={:?} PL={:?}",
            s30.gt.as_ref().map(|g| &g.alleles),
            s30.ad,
            s30.pl
        ),
    );
    assert_eq!(s30.ad.as_deref(), Some(&[42u32, 5][..]));
    assert_eq!(s30.pl.as_deref(), Some(&[84u32, 0, 1738][..]));

    let outcome10 = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &args_at_stand(DEFAULT_STAND_EMIT_CONFIDENCE),
    )
    .expect("call10")
    .expect("Some");
    let call = outcome10
        .genotyped_calls
        .iter()
        .find(|c| call_is_gc(c))
        .expect("stand=10 still genotypes G>C");
    let pl = call.genotype.format.pl_as_i32();
    kv("gt_pl_qual_stand10", format!("PL={pl:?}"));
    assert_eq!(pl, vec![21, 0, 1461]);
    assert_eq!(best_pl_index(&call.genotype.format.pl), 1, "GT 0/1");
    let gls = &call.genotype.genotype_log10_likelihoods;
    let af10 = java_emit_af_decision(gls, DEFAULT_STAND_EMIT_CONFIDENCE).expect("af10");
    let af30 = java_emit_af_decision(gls, DEFAULT_STAND_CALL_CONF).expect("af30");
    kv(
        "qual",
        format!(
            "stand10={:.4} stand30_mono={} stand30_qual={:.4}",
            af10.phred_scaled, af30.site_is_monomorphic, af30.phred_scaled
        ),
    );
    assert!((af10.phred_scaled - 13.63).abs() < 0.02);
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
        try_emit_call_region_variants(region, &outcome10, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    assert!(
        recs10
            .iter()
            .any(|r| rec_is(r, TARGET, TARGET_REF, TARGET_ALT)),
        "same GLs emit at threshold 10"
    );

    let mut keep10 = vec![call.clone()];
    let mut cfg10 = HcGenotypingConfig::strict_java();
    cfg10.stand_emit_confidence = DEFAULT_STAND_EMIT_CONFIDENCE;
    filter_genotyped_calls_for_strict_java_emit(
        &mut keep10,
        &region.reads,
        &outcome10.assembly,
        &cfg10,
    )
    .expect("filter10");
    assert!(keep10.iter().any(call_is_gc), "filter at 10 retains G>C");

    let mut drop30 = vec![call.clone()];
    filter_genotyped_calls_for_strict_java_emit(
        &mut drop30,
        &region.reads,
        &outcome10.assembly,
        &HcGenotypingConfig::strict_java(),
    )
    .expect("filter30");
    assert!(
        !drop30.iter().any(call_is_gc),
        "actual strict-Java filter at 30 drops G>C"
    );
    kv("classification", "EMISSION_PREDICATE_DIVERGENCE");
    kv(
        "production_change",
        "HcGenotypingConfig::strict_java stand_emit_confidence=30",
    );
}

#[test]
fn forensic_6r241_same_pl_fails_java_calling_threshold_only() {
    let event = gatk_haplotypecaller::event_map::VariationEvent::from_alleles(
        "20", TARGET, TARGET_REF, TARGET_ALT,
    );
    let pl = [21, 0, 1461];
    let gl: Vec<f64> = pl.iter().map(|&p| (p as f64) / -10.0).collect();
    let fmt =
        gatk_haplotypecaller::genotyping::emit_genotype_format_fields(&gl, &[35, 3]).expect("fmt");
    assert_eq!(best_pl_index(&fmt.pl), 1);
    let af10 = java_emit_af_decision(&gl, DEFAULT_STAND_EMIT_CONFIDENCE).expect("af10");
    let af30 = java_emit_af_decision(&gl, DEFAULT_STAND_CALL_CONF).expect("af30");
    assert!(af10.passes_emit);
    assert!((af10.phred_scaled - 13.6327).abs() < 0.01);
    assert!(af30.site_is_monomorphic);
    assert!(!af30.passes_emit);
    assert!(java_emit_would_pass(&event, &gl, &fmt, DEFAULT_STAND_EMIT_CONFIDENCE, &[]).unwrap());
    assert!(!java_emit_would_pass(&event, &gl, &fmt, DEFAULT_STAND_CALL_CONF, &[]).unwrap());
}
