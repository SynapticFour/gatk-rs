//! 6R.174: INFO DP is Java `Coverage.evidenceCount`, not FORMAT/`DepthPerSampleHC`.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! Target `2:92305634 G/T`. SOR is out of scope. FORMAT DP stays 2.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r174_info_dp_java_coverage_source -- --nocapture --test-threads=1
//! HOLDOUT_6R174=1 cargo test -p gatk-haplotypecaller --test holdout_6r174_info_dp -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::hc_genotyping_engine::java_alignment_read_overlaps_interval;
use gatk_haplotypecaller::read_assembly_filter::{passes_assembly_read, AssemblyReadFilterConfig};
use gatk_haplotypecaller::variant_site_hc_annotations::coverage_evidence_count;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
    DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use rust_htslib::bam::Record;
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92305500-92305850";
const CLOSED_INTERVAL: &str = "2:92316200-92316580";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_305_634;
const CLOSED: u64 = 92_316_347;
const MERGED_REF: &str = "G";
const MERGED_ALT: &str = "T";
const MARGIN: i32 = DEFAULT_INFORMATIVE_READ_OVERLAP_MARGIN;
const FLAG_REVERSE: u16 = 0x10;
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_DUP: u16 = 0x400;
const FLAG_SUPPLEMENTARY: u16 = 0x800;
const THIRD_READ: &str = "H06JUADXX130110:1:1101:10018:4569 FLAG=145 MAPQ=25 strand=rev";
const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R174\t{key}\t{}", value.as_ref());
}

fn info_i32(info: &[InfoValue], key: &str) -> Option<i32> {
    for v in info {
        if let InfoValue::Integer(k, xs) = v {
            if k == key {
                return xs.first().copied();
            }
        }
    }
    None
}

fn info_f64(info: &[InfoValue], key: &str) -> Option<f64> {
    for v in info {
        if let InfoValue::Float(k, xs) = v {
            if k == key {
                return xs.first().copied();
            }
        }
    }
    None
}

fn member_id(rec: &Record) -> String {
    format!(
        "{} FLAG={} MAPQ={} strand={}",
        String::from_utf8_lossy(rec.qname()),
        rec.flags(),
        rec.mapq(),
        if rec.flags() & FLAG_REVERSE != 0 {
            "rev"
        } else {
            "fwd"
        }
    )
}

fn dump_read(prefix: &str, rec: &Record) {
    eprintln!(
        "{prefix}\tQNAME={} FLAG={} pos={} end={} CIGAR={} MAPQ={} strand={} dup={} sec={} sup={} filtered_asm={}",
        String::from_utf8_lossy(rec.qname()),
        rec.flags(),
        rec.pos() + 1,
        gatk_haplotypecaller::read_unclip::alignment_end_1based(rec),
        rec.cigar(),
        rec.mapq(),
        if rec.flags() & FLAG_REVERSE != 0 {
            "rev"
        } else {
            "fwd"
        },
        rec.flags() & FLAG_DUP != 0,
        rec.flags() & FLAG_SECONDARY != 0,
        rec.flags() & FLAG_SUPPLEMENTARY != 0,
        !passes_assembly_read(rec, &AssemblyReadFilterConfig::gatk_defaults()),
    );
}

#[test]
fn forensic_6r174_java_coverage_contract() {
    kv("java_pin", JAVA_PIN);
    let src = include_str!("../src/variant_site_hc_annotations.rs");
    assert!(
        src.contains("Coverage.annotate")
            && src.contains("evidenceCount")
            && src.contains("coverage_evidence_count"),
        "INFO DP must name Java 4.4.0.0 Coverage.evidenceCount"
    );
    assert!(
        src.contains("6R.174: INFO DP is Java `Coverage.evidenceCount`"),
        "production assignment must be 6R.174 Coverage, not FORMAT DP"
    );
    assert!(
        src.contains("coverage_evidence_count("),
        "INFO DP must call coverage_evidence_count when likelihoods are present"
    );
    assert!(!src.contains("92305634"), "no locus-specific DP patch");
    kv(
        "java_coverage",
        "Coverage.annotate → likelihoods.evidenceCount() (StandardAnnotation INFO DP)",
    );
    kv(
        "java_format_dp",
        "DepthPerSampleHC → informative bestAlleles == sum(AD)",
    );
}

#[test]
fn forensic_6r174_info_dp_java_coverage_source() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
        return;
    }
    kv("java_pin", JAVA_PIN);
    kv("target", "2:92305634 G/T");
    kv("java_info_dp_before", "3");
    kv("rust_info_dp_before", "2");

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
        .expect("genotyped G/T");
    let fmt = &call.genotype.format;
    assert_eq!(fmt.ad_as_i32(), vec![0, 2], "FORMAT AD unchanged");
    assert_eq!(fmt.pl_as_i32(), vec![90, 6, 0], "FORMAT PL unchanged");
    assert_eq!(fmt.dp.as_i32(), 2, "FORMAT DP unchanged (DepthPerSampleHC)");
    assert_eq!(fmt.gq.as_i32(), 6, "FORMAT GQ unchanged");

    let mut retain = BTreeSet::new();
    for cell in &outcome.read_likelihoods {
        let idx = cell.read_index.get();
        let Some(r) = outcome.genotyping_reads.get(idx) else {
            continue;
        };
        if java_alignment_read_overlaps_interval(r, TARGET, TARGET, MARGIN) {
            retain.insert(idx);
        }
    }
    kv("retainEvidence_unique", retain.len().to_string());
    let mut members = Vec::new();
    for idx in &retain {
        let rec = &outcome.genotyping_reads[*idx];
        dump_read("RETAIN", rec);
        members.push(member_id(rec));
    }
    assert_eq!(retain.len(), 3, "Java Coverage evidenceCount = 3");
    assert!(
        members.iter().any(|m| m == THIRD_READ),
        "third INFO DP read must be {THIRD_READ}, got {members:?}"
    );
    kv("third_read", THIRD_READ);
    kv(
        "third_read_role",
        "present in post-filter overlapping sampleEvidence; not counted in FORMAT DP/AD=0,2",
    );

    let counted = coverage_evidence_count(
        &outcome.genotyping_reads,
        &outcome.read_likelihoods,
        TARGET,
        TARGET,
        MARGIN,
    );
    kv("coverage_evidence_count", counted.to_string());
    assert_eq!(counted, 3);
    assert_ne!(
        counted,
        fmt.dp.as_i32(),
        "INFO Coverage count must not collapse to FORMAT DP"
    );

    let mut walker_filtered = 0usize;
    for rec in &covering.reads {
        if !java_alignment_read_overlaps_interval(rec, TARGET, TARGET, MARGIN) {
            continue;
        }
        dump_read("WALKER_OVERLAP", rec);
        if !passes_assembly_read(rec, &AssemblyReadFilterConfig::gatk_defaults()) {
            walker_filtered += 1;
            kv("walker_filtered_overlap", member_id(rec));
        }
    }
    kv(
        "addEvidence_filtered_overlap_candidates",
        walker_filtered.to_string(),
    );

    let emitted =
        try_emit_call_region_variants(covering, &outcome, "NA12878", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let rec = emitted
        .iter()
        .find(|r| r.position == TARGET)
        .expect("emitted");
    let sample = rec.samples.first().expect("sample");
    assert_eq!(sample.dp.map(|d| d as i32), Some(2), "FORMAT DP still 2");
    assert_eq!(
        info_i32(&rec.info, "DP"),
        Some(3),
        "INFO DP = Java Coverage 3"
    );
    kv("rust_info_dp_after", "3");
    kv("format_dp_after", "2");

    // Closed 6R.164/168 target: FORMAT/FS/SOR/MQ unchanged.
    let specs2 = parse_intervals_cli_string(&dict, CLOSED_INTERVAL).expect("closed interval");
    let walk2 = traverse_assembly_region_walker(
        &dict,
        &specs2,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("walk2");
    let regions2 = flatten_assembly_regions(&walk2);
    let covering2 = regions2
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED
                && r.end.get() >= CLOSED
        })
        .expect("closed covering");
    let outcome2 = HaplotypeCallerEngine::call_region(
        covering2,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call2")
    .expect("outcome2");
    let call2 = outcome2
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED)
                && c.event.ref_allele == "G"
                && c.event.alt_allele == "A"
        })
        .expect("closed G/A");
    assert_eq!(call2.genotype.format.ad_as_i32(), vec![0, 3]);
    assert_eq!(call2.genotype.format.dp.as_i32(), 3);
    assert_eq!(call2.genotype.format.pl_as_i32(), vec![135, 9, 0]);
    let emitted2 = try_emit_call_region_variants(
        covering2,
        &outcome2,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("emit2");
    let rec2 = emitted2
        .iter()
        .find(|r| r.position == CLOSED)
        .expect("closed emit");
    assert_eq!(info_i32(&rec2.info, "DP"), Some(3));
    let fs = info_f64(&rec2.info, "FS").unwrap_or(-1.0);
    let sor = info_f64(&rec2.info, "SOR").unwrap_or(-1.0);
    let mq = info_f64(&rec2.info, "MQ").unwrap_or(-1.0);
    assert!((fs - 0.0).abs() < 1e-6, "FS stays 0");
    assert!((sor - 1.179).abs() < 0.002, "SOR stays 1.179, got {sor}");
    assert!((mq - 40.25).abs() < 1e-9, "MQ stays 40.25, got {mq}");
    kv("closed_92316347", "FORMAT/FS/SOR/MQ unchanged");
}
