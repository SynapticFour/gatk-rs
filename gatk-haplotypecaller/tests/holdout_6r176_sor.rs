//! 6R.176 live: SOR at `2:92305634 G/T` is Java 0.693; INFO DP stays 3.
//! Skipped unless `HOLDOUT_6R176=1`.
//!
//! ```text
//! HOLDOUT_6R176=1 cargo test -p gatk-haplotypecaller --test holdout_6r176_sor -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::variant_site_hc_annotations::{
    strand_bias_contingency_table, HcStrandBiasLikelihoods, STRAND_ODDS_RATIO_MIN_COUNT,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, HcGenotypingConfig, ReadFilterParams, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92305500-92305850";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_305_634;
const MERGED_REF: &str = "G";
const MERGED_ALT: &str = "T";

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R176\t{key}\t{}", value.as_ref());
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

#[test]
fn holdout_6r176_sor() {
    if std::env::var("HOLDOUT_6R176").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R176=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
        return;
    }
    kv(
        "java_pin",
        "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77",
    );
    kv("target", "2:92305634 G/T");
    kv(
        "production_change",
        "PairHMM membership ← MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE; addEvidence(0) for Coverage",
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
        .expect("call");
    assert_eq!(call.genotype.format.dp.as_i32(), 2);
    assert_eq!(call.genotype.format.ad_as_i32(), vec![0, 2]);
    assert_eq!(call.genotype.format.pl_as_i32(), vec![90, 6, 0]);
    assert_eq!(call.genotype.format.gq.as_i32(), 6);

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
        .expect("record");
    kv(
        "info_dp",
        info_i32(&rec.info, "DP").unwrap_or(-1).to_string(),
    );
    kv(
        "format_dp",
        rec.samples
            .first()
            .and_then(|s| s.dp)
            .map(|d| d.to_string())
            .unwrap_or_default(),
    );
    let sor = info_f64(&rec.info, "SOR").unwrap_or(0.0);
    kv("sor", format!("{sor}"));
    assert_eq!(info_i32(&rec.info, "DP"), Some(3));
    assert_eq!(
        rec.samples.first().and_then(|s| s.dp).map(|d| d as i32),
        Some(2)
    );
    assert!(
        (sor - 0.6931471805599453).abs() < 1e-6 || (sor - 0.693).abs() < 0.002,
        "SOR must be Java 0.693, got {sor}"
    );
    assert!(
        (rec.quality.unwrap_or(0.0) - 78.32).abs() < 0.02,
        "QUAL stays 78.32"
    );
}
