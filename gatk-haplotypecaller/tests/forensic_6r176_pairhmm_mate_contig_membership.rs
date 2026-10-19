//! 6R.176: PairHMM mate-contig membership at `2:92305634 G/T`.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! FLAG=145 (mate on contig 9) stays in Coverage / INFO DP=3 and is not PairHMM-scored.
//! SOR table `[0,0;1,1]` → 0.693. FORMAT stays 6R.174-closed.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r176_pairhmm_mate_contig_membership -- --nocapture --test-threads=1
//! HOLDOUT_6R176=1 cargo test -p gatk-haplotypecaller --test holdout_6r176_sor -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::variant_site_hc_annotations::{
    coverage_evidence_count, strand_bias_contingency_table, HcStrandBiasLikelihoods,
    STRAND_ODDS_RATIO_MIN_COUNT,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, HcGenotypingConfig, ReadFilterParams, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92305500-92305850";
const CLOSED_INTERVAL: &str = "2:92316200-92316580";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_305_634;
const CLOSED: u64 = 92_316_347;
const MERGED_REF: &str = "G";
const MERGED_ALT: &str = "T";
const THIRD_QNAME: &str = "H06JUADXX130110:1:1101:10018:4569";
const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R176\t{key}\t{}", value.as_ref());
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

fn qname(rec: &rust_htslib::bam::Record) -> String {
    String::from_utf8_lossy(rec.qname()).into_owned()
}

#[test]
fn forensic_6r176_java_mate_contig_pairhmm_contract() {
    kv("java_pin", JAVA_PIN);
    kv(
        "java_path",
        "filterNonPassingReads MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE → PairHMM skip; addEvidence(..., 0) for Coverage",
    );
    let src = include_str!("../src/engine_likelihoods.rs");
    assert!(
        src.contains("score_pairhmm_from_records_java_mate_contig")
            && src.contains("MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE")
            && src.contains("addEvidence"),
        "PairHMM construction must name the Java mate-contig gate"
    );
    assert!(
        !src.contains("92305634") && !src.contains(THIRD_QNAME),
        "no locus- or QNAME-specific PairHMM patch"
    );
    let coverage = include_str!("../src/variant_site_hc_annotations.rs");
    assert!(
        coverage.contains("6R.174: INFO DP is Java `Coverage.evidenceCount`"),
        "INFO DP assignment stays 6R.174 Coverage"
    );
}

#[test]
fn forensic_6r176_pairhmm_mate_contig_membership() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
        return;
    }
    kv("java_pin", JAVA_PIN);
    kv("target", "2:92305634 G/T");
    kv(
        "production_change",
        "compute_region_read_likelihoods mate-contig PairHMM membership; addEvidence(0) cells",
    );

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
        .expect("covering");

    let bam_flag145 = covering
        .reads
        .iter()
        .find(|r| qname(r) == THIRD_QNAME && r.flags() == 145)
        .expect("FLAG=145 in walker evidence");
    kv(
        "flag145_bam",
        format!(
            "FLAG={} MAPQ={} tid={} mtid={} contig={} mate_contig={} start={} CIGAR={}",
            bam_flag145.flags(),
            bam_flag145.mapq(),
            bam_flag145.tid(),
            bam_flag145.mtid(),
            dict.contig_records()
                .get(bam_flag145.tid() as usize)
                .map(|s| s.name.as_str())
                .unwrap_or("?"),
            dict.contig_records()
                .get(bam_flag145.mtid() as usize)
                .map(|s| s.name.as_str())
                .unwrap_or("?"),
            bam_flag145.pos() + 1,
            bam_flag145.cigar()
        ),
    );
    assert_eq!(bam_flag145.tid(), 1, "read is on contig 2 (tid=1)");
    assert_eq!(bam_flag145.mtid(), 8, "mapped mate is on contig 9 (tid=8)");
    assert!(bam_flag145.is_paired(), "FLAG=145 is paired");
    assert!(
        !bam_flag145.is_mate_unmapped(),
        "FLAG=145 mate is mapped (not 0x8)"
    );
    assert_ne!(
        bam_flag145.tid(),
        bam_flag145.mtid(),
        "Java MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE must fail"
    );

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
    assert_eq!(fmt.ad_as_i32(), vec![0, 2]);
    assert_eq!(fmt.pl_as_i32(), vec![90, 6, 0]);
    assert_eq!(fmt.dp.as_i32(), 2);
    assert_eq!(fmt.gq.as_i32(), 6);

    let cfg = HcGenotypingConfig::default();
    let (full_ref, full_pad) = outcome.assembly.event_map_reference();
    let apply_pad = outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.as_ref().map(|g| g.start.get()))
        .unwrap_or(full_pad);
    let evidence = HcStrandBiasLikelihoods {
        reads: &outcome.genotyping_reads,
        likelihoods: &outcome.read_likelihoods,
        haplotypes: &outcome.assembly.haplotypes,
        contig: &covering.contig,
        ref_bytes: outcome.assembly.reference_bases(),
        pad_start_1based: apply_pad,
        full_ref_bytes: full_ref,
        full_pad_1based: full_pad,
        max_mnp_distance: outcome.assembly.max_mnp_distance(),
        emit_spanning_dels: !cfg.disable_spanning_event_genotyping,
    };

    let coverage_only =
        gatk_haplotypecaller::with_coverage_only_failed_mates(|mates| mates.to_vec());
    let third_reads: Vec<_> = coverage_only
        .iter()
        .filter(|r| qname(r) == THIRD_QNAME && r.flags() == 145)
        .collect();
    assert_eq!(
        third_reads.len(),
        1,
        "FLAG=145 remains in the coverage-only list"
    );
    let rec = third_reads[0];
    assert!(
        outcome
            .genotyping_reads
            .iter()
            .all(|g| qname(g) != THIRD_QNAME || g.flags() != 145),
        "FLAG=145 must stay out of genotyping evidence"
    );
    kv(
        "flag145",
        format!(
            "QNAME={} FLAG={} MAPQ={} tid={} mtid={} start={} CIGAR={}",
            qname(rec),
            rec.flags(),
            rec.mapq(),
            rec.tid(),
            rec.mtid(),
            rec.pos() + 1,
            rec.cigar()
        ),
    );

    let third_cells: Vec<f64> = evidence
        .likelihoods
        .iter()
        .filter(|c| {
            evidence
                .reads
                .get(c.read_index.get())
                .is_some_and(|r| qname(r) == THIRD_QNAME && r.flags() == 145)
        })
        .map(|c| c.log10_likelihood)
        .collect();
    assert!(
        third_cells.is_empty(),
        "FLAG=145 must not have a PairHMM or ordinary likelihood row, got {third_cells:?}"
    );
    kv("flag145_pairhmm", "no_likelihood_row");

    let counted = coverage_evidence_count(
        evidence.reads,
        evidence.likelihoods,
        TARGET,
        TARGET,
        cfg.informative_read_overlap_margin,
    );
    kv("coverage_evidence_count", counted.to_string());
    assert_eq!(counted, 3, "INFO DP Coverage still includes FLAG=145");

    let table = strand_bias_contingency_table(
        &evidence,
        TARGET,
        MERGED_REF,
        MERGED_ALT,
        STRAND_ODDS_RATIO_MIN_COUNT,
    );
    kv(
        "sor_table",
        format!("[{},{};{},{}]", table.0, table.1, table.2, table.3),
    );
    assert_eq!(table, (0, 0, 1, 1));

    let emitted =
        try_emit_call_region_variants(covering, &outcome, "NA12878", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let rec = emitted
        .iter()
        .find(|r| r.position == TARGET)
        .expect("emitted");
    let sample = rec.samples.first().expect("sample");
    assert_eq!(info_i32(&rec.info, "DP"), Some(3));
    assert_eq!(sample.dp.map(|d| d as i32), Some(2));
    let sor = info_f64(&rec.info, "SOR").unwrap_or(-1.0);
    kv("sor", format!("{sor}"));
    assert!(
        (sor - 0.6931471805599453).abs() < 1e-6 || (sor - 0.693).abs() < 0.002,
        "SOR must be Java 0.693, got {sor}"
    );
    assert!(
        (rec.quality.unwrap_or(0.0) - 78.32).abs() < 0.02,
        "QUAL stays 78.32"
    );
    kv(
        "remaining_info",
        format!(
            "MQ={:?} FS={:?} ReadPosRankSum={:?} InbreedingCoeff={:?}",
            info_f64(&rec.info, "MQ"),
            info_f64(&rec.info, "FS"),
            info_f64(&rec.info, "ReadPosRankSum"),
            info_f64(&rec.info, "InbreedingCoeff")
        ),
    );

    let specs2 = parse_intervals_cli_string(&dict, CLOSED_INTERVAL).expect("closed");
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
    let sor2 = info_f64(&rec2.info, "SOR").unwrap_or(-1.0);
    let mq2 = info_f64(&rec2.info, "MQ").unwrap_or(-1.0);
    assert!((sor2 - 1.179).abs() < 0.002, "6R.166 SOR stays 1.179");
    assert!((mq2 - 40.25).abs() < 1e-9, "6R.168 MQ stays 40.25");
    kv("closed_92316347", "FS/SOR/MQ unchanged");
}
