//! 6R.167: MQ consumes Java `RMSMappingQuality.calculateRawData` sampleEvidence.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`, GKL 0.8.8.
//! Target `2:92316347 G/A`. FS/SOR stay 6R.166-closed. FORMAT/QUAL stay 6R.164-closed.
//! Aggregation remains the pre-existing arithmetic mean; Java RMS is 6R.168.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r167_mq_java_evidence_source -- --test-threads=1 --nocapture
//! HOLDOUT_6R167=1 cargo test -p gatk-haplotypecaller --test holdout_6r167_mq_evidence -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::fragment_overlap::read_base_at_ref_coord_1based;
use gatk_haplotypecaller::read_model::MAPPING_QUALITY_UNAVAILABLE;
use gatk_haplotypecaller::variant_site_hc_annotations::{
    rms_mapping_quality_sample_mapqs, rms_mapping_quality_sample_reads,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use rust_htslib::bam::Record;
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92316200-92316580";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_316_347;
const MERGED_REF: &str = "G";
const MERGED_ALT: &str = "A";
const FLAG_REVERSE: u16 = 0x10;
const FLAG_SECONDARY: u16 = 0x100;
const FLAG_DUP: u16 = 0x400;
const FLAG_SUPPLEMENTARY: u16 = 0x800;
/// Java `RMSMappingQuality.calculateRawData` members (QNAME + FLAG + MAPQ).
const JAVA_MQ_MEMBERS: &[&str] = &[
    "H06HDADXX130110:1:1101:10034:45116 FLAG=99 MAPQ=47 strand=fwd",
    "H06HDADXX130110:2:1101:10025:49248 FLAG=99 MAPQ=47 strand=fwd",
    "H06HDADXX130110:2:1101:10046:78083 FLAG=147 MAPQ=21 strand=rev",
];
/// 6R.165 ALT-pileup extras that are not in post-filter `sampleEvidence`.
const RUST_PILEUP_ONLY: &[&str] = &[
    "H06HDADXX130110:2:1101:10025:49248 FLAG=147 MAPQ=46 strand=rev",
    "H06HDADXX130110:2:1101:10046:78083 FLAG=99 MAPQ=21 strand=fwd",
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

fn member_id(rec: &Record) -> String {
    let reverse = rec.flags() & FLAG_REVERSE != 0;
    format!(
        "{} FLAG={} MAPQ={} strand={}",
        String::from_utf8_lossy(rec.qname()),
        rec.flags(),
        rec.mapq(),
        if reverse { "rev" } else { "fwd" }
    )
}

fn dump_read(prefix: &str, rec: &Record) {
    eprintln!(
        "{prefix}\tQNAME={} FLAG={} pos={} CIGAR={} MAPQ={} strand={} dup={} sec={} sup={}",
        String::from_utf8_lossy(rec.qname()),
        rec.flags(),
        rec.pos() + 1,
        rec.cigar(),
        rec.mapq(),
        if rec.flags() & FLAG_REVERSE != 0 {
            "rev"
        } else {
            "fwd"
        },
        rec.flags() & FLAG_DUP != 0,
        rec.flags() & FLAG_SECONDARY != 0,
        rec.flags() & FLAG_SUPPLEMENTARY != 0
    );
}

fn java_rms(mqs: &[u8]) -> (usize, u64, f64) {
    let n = mqs.len();
    let sum_sq: u64 = mqs.iter().map(|&m| u64::from(m) * u64::from(m)).sum();
    let rms = if n == 0 {
        0.0
    } else {
        (sum_sq as f64 / n as f64).sqrt()
    };
    (n, sum_sq, rms)
}

fn mean(mqs: &[u8]) -> f64 {
    if mqs.is_empty() {
        0.0
    } else {
        mqs.iter().map(|&m| u64::from(m)).sum::<u64>() as f64 / mqs.len() as f64
    }
}

#[test]
fn forensic_6r167_java_rmsmappingquality_contract() {
    assert_eq!(MAPPING_QUALITY_UNAVAILABLE, 255);
    let src = include_str!("../src/variant_site_hc_annotations.rs");
    assert!(
        src.contains("RMSMappingQuality.calculateRawData")
            && src.contains("sampleEvidence")
            && src.contains("MAPPING_QUALITY_UNAVAILABLE"),
        "production MQ must name Java 4.4.0.0 calculateRawData / sampleEvidence"
    );
    assert!(
        src.contains("rms_mapping_quality_sample_mapqs")
            && src.contains("6R.167: MQ membership is Java"),
        "6R.167 sampleEvidence membership must stay"
    );
    assert!(
        src.contains("mq_rms_of_sample_evidence") && src.contains("makeFinalizedAnnotationString"),
        "6R.168 owns the RMS formula"
    );
    assert!(!src.contains("92316347"), "no locus-specific MQ patch");
}

#[test]
fn forensic_6r167_mq_java_evidence_source() {
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

    let alt_b = MERGED_ALT.as_bytes()[0];
    let mut pileup_alt = Vec::new();
    for rec in covering.reads.iter().map(|r| r.as_ref()) {
        let Some(base) = read_base_at_ref_coord_1based(rec, TARGET as i32) else {
            continue;
        };
        if !base.eq_ignore_ascii_case(&alt_b) {
            continue;
        }
        dump_read("6R167\trust_before_alt_pileup", rec);
        pileup_alt.push(member_id(rec));
    }
    let pileup_set: BTreeSet<_> = pileup_alt.iter().cloned().collect();
    let java_set: BTreeSet<_> = JAVA_MQ_MEMBERS.iter().map(|s| (*s).to_string()).collect();
    let pileup_only: BTreeSet<_> = pileup_set.difference(&java_set).cloned().collect();
    eprintln!("6R167\trust_before_mapqs\t{:?}", {
        covering
            .reads
            .iter()
            .filter_map(|r| {
                read_base_at_ref_coord_1based(r, TARGET as i32)
                    .filter(|b| b.eq_ignore_ascii_case(&alt_b))
                    .map(|_| r.mapq())
            })
            .collect::<Vec<_>>()
    });
    assert_eq!(
        pileup_alt.len(),
        5,
        "6R.165 ALT pileup must still be n=5 before the evidence switch"
    );
    for extra in RUST_PILEUP_ONLY {
        assert!(
            pileup_only.iter().any(|m| m == extra),
            "Rust-only pileup member missing: {extra}; have {pileup_only:?}"
        );
    }

    let sample_reads =
        rms_mapping_quality_sample_reads(&outcome.genotyping_reads, &outcome.read_likelihoods);
    for rec in &sample_reads {
        dump_read("6R167\tjava_sample_evidence", rec.as_ref());
        assert_ne!(
            rec.mapq(),
            MAPPING_QUALITY_UNAVAILABLE,
            "this site has no MQ=255 evidence"
        );
        assert_eq!(
            rec.flags() & (FLAG_DUP | FLAG_SECONDARY | FLAG_SUPPLEMENTARY),
            0
        );
    }
    let after_ids: Vec<String> = sample_reads.iter().map(|r| member_id(r)).collect();
    let after_set: BTreeSet<_> = after_ids.iter().cloned().collect();
    assert_eq!(
        after_set, java_set,
        "Rust sampleEvidence membership must match Java RMSMappingQuality"
    );
    assert!(
        outcome.genotyping_reads.len() > after_ids.len(),
        "genotyping_reads is a superset; MQ must not consume the full overlap list"
    );
    let after_mapqs =
        rms_mapping_quality_sample_mapqs(&outcome.genotyping_reads, &outcome.read_likelihoods);
    eprintln!("6R167\trust_after_mapqs\t{after_mapqs:?}");
    let mut sorted_after = after_mapqs.clone();
    sorted_after.sort_unstable();
    assert_eq!(sorted_after, vec![21, 47, 47]);

    let (n, sum_sq, rms) = java_rms(&after_mapqs);
    let rust_mean = mean(&after_mapqs);
    eprintln!("6R167\tjava_calc\tn={n} sum_sq={sum_sq} rms={rms:.5} printed={rms:.2}");
    eprintln!("6R167\trust_mean\t{rust_mean}");
    assert_eq!(n, 3);
    assert_eq!(sum_sq, 4859);
    assert!((rms - (4859.0_f64 / 3.0).sqrt()).abs() < 1e-12);
    assert!((rms - 40.25).abs() < 0.005);
    assert!((rust_mean - 115.0 / 3.0).abs() < 1e-12);
    assert!(
        (rust_mean - 40.25).abs() > 1.0,
        "membership match must not hide the remaining mean-vs-RMS split"
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
        "6R167\temit\tGT={gt} AD={ad} PL={pl} QUAL={:?} FS={fs} SOR={sor} MQ={mq}",
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
        "after 6R.168, production MQ is Java 40.25, got {mq}"
    );
    assert!(
        (mq - 36.4).abs() > 1.0,
        "production MQ must leave the 6R.165 ALT-pileup mean"
    );
}
