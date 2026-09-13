//! 6R.166: FS/SOR consume Java `getContingencyTable` evidence, not `region.reads` pileup.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`, GKL 0.8.8.
//! Target `2:92316347 G/A`. FORMAT/QUAL stay 6R.164-closed.
//! MQ evidence moved in 6R.167; this file still pins FS/SOR.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r166_fs_sor_java_evidence_source -- --test-threads=1 --nocapture
//! HOLDOUT_6R166=1 cargo test -p gatk-haplotypecaller --test holdout_6r166_info_annotation -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{build_per_haplotype_variation_events, VariationEvent};
use gatk_haplotypecaller::hc_allele_mapping::create_allele_mapper_with_events;
use gatk_haplotypecaller::read_realignment::LOG_10_INFORMATIVE_THRESHOLD;
use gatk_haplotypecaller::variant_site_hc_annotations::{
    strand_bias_contingency_table, HcStrandBiasLikelihoods, FISHER_STRAND_MIN_COUNT,
    STRAND_ODDS_RATIO_MIN_COUNT,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, marginalize_rows_to_biallelic_alleles,
    region_likelihoods_to_rows, traverse_assembly_region_walker, try_emit_call_region_variants,
    AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition, HaplotypeCallerEngine,
    HcGenotypingConfig, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92316200-92316580";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_316_347;
const MERGED_REF: &str = "G";
const MERGED_ALT: &str = "A";
const FLAG_REVERSE: u16 = 0x10;
/// Java 4.4.0.0 informative best-allele members at this site (post-filter AlleleLikelihoods).
const JAVA_INFORMATIVE_MEMBERS: &[&str] = &[
    "H06HDADXX130110:1:1101:10034:45116 FLAG=99 strand=fwd allele=ALT",
    "H06HDADXX130110:2:1101:10025:49248 FLAG=99 strand=fwd allele=ALT",
    "H06HDADXX130110:2:1101:10046:78083 FLAG=147 strand=rev allele=ALT",
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

#[test]
fn forensic_6r166_java_get_contingency_table_contract() {
    assert_eq!(FISHER_STRAND_MIN_COUNT, 2);
    assert_eq!(STRAND_ODDS_RATIO_MIN_COUNT, 0);
    let src = include_str!("../src/variant_site_hc_annotations.rs");
    assert!(
        src.contains("StrandBiasTest.getContingencyTable")
            && src.contains("bestAllelesBreakingTies")
            && src.contains("isInformative"),
        "production FS/SOR table must name the Java 4.4.0.0 methods"
    );
    assert!(
        src.contains("6R.166: FS/SOR use Java `getContingencyTable`")
            && src.contains("6R.167: MQ membership is Java"),
        "FS/SOR stay on getContingencyTable; MQ membership stays 6R.167"
    );
    assert!(!src.contains("92316347"), "no locus-specific FS/SOR patch");
}

#[test]
fn forensic_6r166_fs_sor_java_evidence_source() {
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

    let (full_ref, full_pad) = outcome.assembly.event_map_reference();
    let apply_pad = outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.as_ref().map(|g| g.start.get()))
        .unwrap_or(full_pad);
    let cfg = HcGenotypingConfig::default();
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
    let hap_cache = build_per_haplotype_variation_events(
        evidence.haplotypes,
        evidence.full_ref_bytes,
        evidence.full_pad_1based,
        evidence.max_mnp_distance,
        evidence.contig,
    );
    let mapping = create_allele_mapper_with_events(
        &VariationEvent::from_alleles(evidence.contig, TARGET, MERGED_REF, MERGED_ALT),
        TARGET,
        evidence.haplotypes,
        evidence.pad_start_1based,
        evidence.ref_bytes,
        evidence.max_mnp_distance,
        evidence.emit_spanning_dels,
        Some(&hap_cache),
    );
    let rows = region_likelihoods_to_rows(evidence.likelihoods, evidence.haplotypes.len());
    let marg = marginalize_rows_to_biallelic_alleles(
        &rows,
        &mapping.ref_haplotype_indices,
        &mapping.alt_haplotype_indices,
    );
    let mut members = Vec::new();
    for row in &marg {
        let Some(rec) = evidence.reads.get(row.read_index) else {
            continue;
        };
        let lls = &row.haplotype_log10_likelihoods;
        let ll_ref = lls.first().copied().unwrap_or(f64::NEG_INFINITY);
        let ll_alt = lls.get(1).copied().unwrap_or(f64::NEG_INFINITY);
        if !ll_ref.is_finite() && !ll_alt.is_finite() {
            continue;
        }
        let (best_is_ref, best, second) = if ll_ref > ll_alt {
            (true, ll_ref, ll_alt)
        } else if ll_alt > ll_ref {
            (false, ll_alt, ll_ref)
        } else {
            (true, ll_ref, ll_alt)
        };
        let gap = if second.is_finite() {
            best - second
        } else {
            f64::INFINITY
        };
        if gap <= LOG_10_INFORMATIVE_THRESHOLD {
            continue;
        }
        let reverse = rec.flags() & FLAG_REVERSE != 0;
        let member = format!(
            "{} FLAG={} strand={} allele={}",
            String::from_utf8_lossy(rec.qname()),
            rec.flags(),
            if reverse { "rev" } else { "fwd" },
            if best_is_ref { "REF" } else { "ALT" }
        );
        eprintln!("6R166\tinformative_read\t{member} gap={gap:.4}");
        members.push(member);
    }
    assert_eq!(
        members, JAVA_INFORMATIVE_MEMBERS,
        "Java-equivalent informative membership, strand, and REF/ALT classification"
    );
    assert!(
        members.iter().all(|m| m.contains("allele=ALT")),
        "this site has no informative REF best allele"
    );

    let rust_fs = strand_bias_contingency_table(
        &evidence,
        TARGET,
        MERGED_REF,
        MERGED_ALT,
        FISHER_STRAND_MIN_COUNT,
    );
    let rust_sor = strand_bias_contingency_table(
        &evidence,
        TARGET,
        MERGED_REF,
        MERGED_ALT,
        STRAND_ODDS_RATIO_MIN_COUNT,
    );
    eprintln!("6R166\trust_table\tFS={:?} SOR={:?}", rust_fs, rust_sor);
    assert_eq!(rust_fs, (0, 0, 2, 1), "Java-equivalent FS table [0,0;2,1]");
    assert_eq!(
        rust_sor,
        (0, 0, 2, 1),
        "at this site SOR MIN_COUNT=0 admits the same three informative reads"
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
        "6R166\temit\tGT={gt} AD={ad} PL={pl} QUAL={:?} FS={fs} SOR={sor} MQ={mq}",
        rec.quality
    );
    assert_eq!(gt, "1/1");
    assert_eq!(ad, "0,3");
    assert_eq!(pl, "135,9,0");
    assert!(
        (rec.quality.unwrap_or(0.0) - 121.84).abs() < 0.02,
        "QUAL must stay 121.84, got {:?}",
        rec.quality
    );
    assert!(fs < 0.02, "FS must be Java 0.000, got {fs}");
    assert!(
        (sor - 1.1786549963416462).abs() < 1e-9 || (sor - 1.179).abs() < 0.001,
        "SOR must be Java 1.17865… (prints 1.179), got {sor}"
    );
    assert!(
        (mq - 40.25).abs() < 1e-12,
        "after 6R.168, MQ is Java RMS 40.25, got {mq}"
    );
}
