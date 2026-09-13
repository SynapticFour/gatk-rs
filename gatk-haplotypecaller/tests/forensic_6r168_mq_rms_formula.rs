//! 6R.168: MQ aggregation is Java `makeFinalizedAnnotationString` RMS, not a mean.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`, GKL 0.8.8.
//! Target `2:92316347 G/A`. 6R.167 sampleEvidence membership stays closed.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r168_mq_rms_formula -- --test-threads=1 --nocapture
//! HOLDOUT_6R168=1 cargo test -p gatk-haplotypecaller --test holdout_6r168_mq_rms -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::read_model::MAPPING_QUALITY_UNAVAILABLE;
use gatk_haplotypecaller::variant_site_hc_annotations::{
    rms_mapping_quality_raw, rms_mapping_quality_sample_mapqs, rms_mapping_quality_sample_reads,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92316200-92316580";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_316_347;
const MERGED_REF: &str = "G";
const MERGED_ALT: &str = "A";
const FLAG_REVERSE: u16 = 0x10;
const JAVA_MQ_MEMBERS: &[&str] = &[
    "H06HDADXX130110:1:1101:10034:45116 FLAG=99 MAPQ=47 strand=fwd",
    "H06HDADXX130110:2:1101:10025:49248 FLAG=99 MAPQ=47 strand=fwd",
    "H06HDADXX130110:2:1101:10046:78083 FLAG=147 MAPQ=21 strand=rev",
];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn info_f64(info: &[InfoValue], key: &str) -> f64 {
    for v in info {
        if let InfoValue::Float(k, xs) = v {
            if k == key {
                return xs.first().copied().unwrap_or(0.0);
            }
        }
    }
    0.0
}

fn member_id(rec: &rust_htslib::bam::Record) -> String {
    let reverse = rec.flags() & FLAG_REVERSE != 0;
    format!(
        "{} FLAG={} MAPQ={} strand={}",
        String::from_utf8_lossy(rec.qname()),
        rec.flags(),
        rec.mapq(),
        if reverse { "rev" } else { "fwd" }
    )
}

#[test]
fn forensic_6r168_java_rms_formula_contract() {
    assert_eq!(MAPPING_QUALITY_UNAVAILABLE, 255);
    let src = include_str!("../src/variant_site_hc_annotations.rs");
    assert!(
        src.contains("makeFinalizedAnnotationString")
            && src.contains("mq_rms_of_sample_evidence")
            && src.contains("java_finalize_mq"),
        "production MQ must use Java RMS + %.2f, not a mean"
    );
    assert!(
        !src.contains("mq_mean_of_sample_evidence"),
        "arithmetic-mean MQ aggregation must be gone"
    );
    assert!(
        src.contains("rms_mapping_quality_sample_mapqs")
            && src.contains("6R.167: MQ membership is Java"),
        "6R.167 sampleEvidence membership must stay"
    );
    assert!(!src.contains("92316347"), "no locus-specific MQ patch");
    let (n, sum_sq, rms) = rms_mapping_quality_raw(&[47, 47, 21]).expect("raw");
    assert_eq!(n, 3);
    assert_eq!(sum_sq, 4859);
    assert!((rms - (4859.0_f64 / 3.0).sqrt()).abs() < 1e-12);
}

#[test]
fn forensic_6r168_mq_rms_formula() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
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
    let covering = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("ActiveFull covering target");

    let outcome = HaplotypeCallerEngine::call_region(
        covering,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("outcome");

    let call = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(TARGET)
                && c.event.ref_allele == MERGED_REF
                && c.event.alt_allele == MERGED_ALT
        })
        .expect("genotyped G/A");
    let fmt = &call.genotype.format;
    assert_eq!(fmt.ad_as_i32(), vec![0, 3]);
    assert_eq!(fmt.pl_as_i32(), vec![135, 9, 0]);
    assert_eq!(fmt.dp.as_i32(), 3);
    assert_eq!(fmt.gq.as_i32(), 9);

    let members: BTreeSet<String> =
        rms_mapping_quality_sample_reads(&outcome.genotyping_reads, &outcome.read_likelihoods)
            .into_iter()
            .map(|r| member_id(r))
            .collect();
    let java_set: BTreeSet<_> = JAVA_MQ_MEMBERS.iter().map(|s| (*s).to_string()).collect();
    assert_eq!(
        members, java_set,
        "6R.167 sampleEvidence membership must be unchanged"
    );

    let mapqs =
        rms_mapping_quality_sample_mapqs(&outcome.genotyping_reads, &outcome.read_likelihoods);
    let mut sorted = mapqs.clone();
    sorted.sort_unstable();
    assert_eq!(sorted, vec![21, 47, 47]);
    eprintln!("6R168\tinput_mapqs\t{mapqs:?}");

    let (n, sum_sq, raw_rms) = rms_mapping_quality_raw(&mapqs).expect("raw RMS");
    eprintln!("6R168\traw\tn={n} sum_sq={sum_sq} rms={raw_rms:.10}");
    assert_eq!(n, 3);
    assert_eq!(sum_sq, 4859);
    assert!((raw_rms - (4859.0_f64 / 3.0).sqrt()).abs() < 1e-12);
    assert!(
        (raw_rms - 38.333333333333336).abs() > 1.0,
        "raw value must be RMS, not the 6R.167 mean"
    );

    let emitted =
        try_emit_call_region_variants(covering, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let rec = emitted
        .iter()
        .find(|r| {
            r.position == TARGET
                && r.reference == MERGED_REF
                && r.alternate.first().map(String::as_str) == Some(MERGED_ALT)
        })
        .expect("emitted G/A");
    let fs = info_f64(&rec.info, "FS");
    let sor = info_f64(&rec.info, "SOR");
    let mq = info_f64(&rec.info, "MQ");
    let sample = rec.samples.first().expect("sample");
    let gt = sample
        .gt
        .as_ref()
        .map(|g| {
            g.alleles
                .iter()
                .map(|a| a.to_string())
                .collect::<Vec<_>>()
                .join("/")
        })
        .unwrap_or_default();
    let ad = sample
        .ad
        .as_ref()
        .map(|v| {
            v.iter()
                .map(|x| x.to_string())
                .collect::<Vec<_>>()
                .join(",")
        })
        .unwrap_or_default();
    let pl = sample
        .pl
        .as_ref()
        .map(|v| {
            v.iter()
                .map(|x| x.to_string())
                .collect::<Vec<_>>()
                .join(",")
        })
        .unwrap_or_default();
    eprintln!(
        "6R168\temit\tGT={gt} AD={ad} PL={pl} QUAL={:?} FS={fs} SOR={sor} MQ={mq}",
        rec.quality
    );
    assert_eq!(gt, "1/1");
    assert_eq!(ad, "0,3");
    assert_eq!(pl, "135,9,0");
    assert!((rec.quality.unwrap_or(0.0) - 121.84).abs() < 0.02);
    assert!(fs < 0.02, "FS must stay Java 0, got {fs}");
    assert!(
        (sor - 1.1786549963416462).abs() < 1e-9 || (sor - 1.179).abs() < 0.001,
        "SOR must stay Java 1.179, got {sor}"
    );
    assert!(
        (mq - 40.25).abs() < 1e-12,
        "production MQ must be Java 40.25, got {mq}"
    );
}
