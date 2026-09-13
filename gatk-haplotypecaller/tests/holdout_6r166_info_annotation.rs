//! 6R.166 live: FS/SOR match Java at `2:92316347 G/A`; FORMAT unchanged.
//! MQ membership moved in 6R.167. Skipped unless `HOLDOUT_6R166=1`.
//!
//! ```text
//! HOLDOUT_6R166=1 cargo test -p gatk-haplotypecaller --test holdout_6r166_info_annotation -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::variant_site_hc_annotations::{
    strand_bias_contingency_table, HcStrandBiasLikelihoods, FISHER_STRAND_MIN_COUNT,
    STRAND_ODDS_RATIO_MIN_COUNT,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, HcGenotypingConfig, ReadFilterParams, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92316200-92316580";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_316_347;
const MERGED_REF: &str = "G";
const MERGED_ALT: &str = "A";

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R166\t{key}\t{}", value.as_ref());
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
fn holdout_6r166_info_annotation() {
    if std::env::var("HOLDOUT_6R166").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R166=1");
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
    kv("variant", format!("2:{TARGET} G/A"));
    kv("production_change", "FS/SOR getContingencyTable evidence");

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
    kv(
        "format",
        format!(
            "AD={:?} PL={:?} DP={} GQ={}",
            call.genotype.format.ad_as_i32(),
            call.genotype.format.pl_as_i32(),
            call.genotype.format.dp.as_i32(),
            call.genotype.format.gq.as_i32()
        ),
    );
    assert_eq!(call.genotype.format.ad_as_i32(), vec![0, 3]);
    assert_eq!(call.genotype.format.pl_as_i32(), vec![135, 9, 0]);

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
    let table = strand_bias_contingency_table(
        &evidence,
        TARGET,
        MERGED_REF,
        MERGED_ALT,
        FISHER_STRAND_MIN_COUNT,
    );
    kv(
        "rust_table",
        format!("[{},{};{},{}]", table.0, table.1, table.2, table.3),
    );
    kv("java_table", "[0,0;2,1]");
    kv("rust_before", "[0,2;3,2]");
    assert_eq!(table, (0, 0, 2, 1));
    assert_eq!(
        strand_bias_contingency_table(
            &evidence,
            TARGET,
            MERGED_REF,
            MERGED_ALT,
            STRAND_ODDS_RATIO_MIN_COUNT,
        ),
        (0, 0, 2, 1)
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
    kv(
        "rust_final",
        format!(
            "FS={:.5} SOR={:.5} MQ={} QUAL={:?}",
            fs, sor, mq, rec.quality
        ),
    );
    kv("java_final", "FS=0 SOR=1.179 MQ=40.25 QUAL=121.84");
    assert!(fs < 0.02, "FS={fs}");
    assert!((sor - 1.179).abs() < 0.001, "SOR={sor}");
    assert!(
        (mq - 40.25).abs() < 1e-12,
        "after 6R.168 MQ is Java RMS 40.25, got {mq}"
    );
    assert!((rec.quality.unwrap_or(0.0) - 121.84).abs() < 0.02);
}
