//! 6R.188 live: `2:92305635 A/G` annotation is the stored-hap loc-loop object (6R.189).
//! FORMAT subset stays AD=0,2. Skipped unless `HOLDOUT_6R188=1`.
//!
//! ```text
//! HOLDOUT_6R188=1 cargo test -p gatk-haplotypecaller --test holdout_6r188_wrong_annotation_source -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92305500-92305850";
const TG_INTERVAL: &str = "2:92307200-92307550";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_305_635;
const CLOSED_SNP: u64 = 92_305_634;
const CLOSED_TG: u64 = 92_307_333;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R188\t{key}\t{}", value.as_ref());
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
fn holdout_6r188_wrong_annotation_source() {
    if std::env::var("HOLDOUT_6R188").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R188=1");
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
    kv(
        "production_change",
        "6R.189 SiteScore stored-hap loc-loop attach",
    );
    kv("target", "2:92305635 A/G");
    kv(
        "classification",
        "A2 closed by 6R.189 — FORMAT subset unchanged; annotation is stored-hap loc-loop n=3",
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

    let stored: BTreeSet<usize> = outcome
        .read_likelihoods
        .iter()
        .map(|c| c.read_index.get())
        .collect();
    assert_eq!(stored.len(), 3);

    let closed = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_SNP)
                && c.event.ref_allele == "G"
                && c.event.alt_allele == "T"
        })
        .expect("closed");
    let target = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(TARGET)
                && c.event.ref_allele == "A"
                && c.event.alt_allele == "G"
        })
        .expect("target");
    assert_eq!(
        closed
            .annotation_likelihoods
            .iter()
            .map(|c| c.read_index.get())
            .collect::<BTreeSet<_>>()
            .len(),
        3,
        "neighbor G/T two-read arm attaches retainEvidence n=3 (6R.205)"
    );
    let attached: BTreeSet<usize> = target
        .annotation_likelihoods
        .iter()
        .map(|c| c.read_index.get())
        .collect();
    assert_eq!(attached.len(), 3, "6R.189 loc-loop n=3");
    assert_eq!(target.genotype.format.ad_as_i32(), vec![0, 2]);
    assert_eq!(target.genotype.format.pl_as_i32(), vec![90, 6, 0]);

    let emitted =
        try_emit_call_region_variants(covering, &outcome, "NA12878", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let closed_rec = emitted
        .iter()
        .find(|r| r.position == CLOSED_SNP && r.reference == "G")
        .expect("closed emit");
    let target_rec = emitted
        .iter()
        .find(|r| r.position == TARGET && r.reference == "A")
        .expect("target emit");
    assert_eq!(info_i32(&closed_rec.info, "DP"), Some(3));
    assert_eq!(info_i32(&target_rec.info, "DP"), Some(3));
    let closed_mq = info_f64(&closed_rec.info, "MQ").unwrap_or(-1.0);
    let target_mq = info_f64(&target_rec.info, "MQ").unwrap_or(-1.0);
    assert!((closed_mq - 41.96).abs() < 0.005);
    assert!((target_mq - 41.96).abs() < 0.005);
    kv(
        "6r186_selected",
        "false; cluster-TG pin is 92307333, not 92305635",
    );

    let tg_specs = parse_intervals_cli_string(&dict, TG_INTERVAL).expect("tg");
    let tg_walk = traverse_assembly_region_walker(
        &dict,
        &tg_specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("tg walk");
    let tg_regions = flatten_assembly_regions(&tg_walk);
    let tg_covering = tg_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_TG
                && r.end.get() >= CLOSED_TG
        })
        .expect("tg covering");
    let tg_outcome = HaplotypeCallerEngine::call_region(
        tg_covering,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("tg call")
    .expect("tg outcome");
    let tg_call = tg_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_TG)
                && c.event.ref_allele == "T"
                && c.event.alt_allele == "G"
        })
        .expect("T/G");
    let tg_ann: BTreeSet<usize> = tg_call
        .annotation_likelihoods
        .iter()
        .map(|c| c.read_index.get())
        .collect();
    assert_eq!(tg_ann.len(), 1, "6R.186 n=1 stays");
}
