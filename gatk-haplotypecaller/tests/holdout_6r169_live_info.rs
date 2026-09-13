//! 6R.169 live: INFO inventory at `2:92316347 G/A`; FORMAT/FS/SOR/MQ unchanged.
//! Skipped unless `HOLDOUT_6R169=1`. Production change: NONE.
//!
//! ```text
//! HOLDOUT_6R169=1 cargo test -p gatk-haplotypecaller --test holdout_6r169_live_info -- --nocapture --test-threads=1
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
    println!("6R169\t{key}\t{}", value.as_ref());
}

fn info_key(v: &InfoValue) -> &str {
    match v {
        InfoValue::Flag(k)
        | InfoValue::Integer(k, _)
        | InfoValue::Float(k, _)
        | InfoValue::String(k, _)
        | InfoValue::Character(k, _) => k.as_str(),
    }
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
fn holdout_6r169_live_info() {
    if std::env::var("HOLDOUT_6R169").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R169=1");
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
    kv("production_change", "NONE");
    kv(
        "note",
        "HOLDOUT_6R160 rust.vcf FORMAT pin is historical (pre-6R.164 AD=0,1); not live INFO",
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
    assert_eq!(call.genotype.format.ad_as_i32(), vec![0, 3]);
    assert_eq!(call.genotype.format.pl_as_i32(), vec![135, 9, 0]);

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

    let keys: BTreeSet<_> = rec.info.iter().map(info_key).collect();
    kv("rust_keys", format!("{keys:?}"));
    kv(
        "java_record_keys",
        "AC,AF,AN,DP,ExcessHet,FS,MLEAC,MLEAF,MQ,QD,SOR",
    );
    kv("rust_only", "ReadPosRankSum");
    let fs = info_f64(&rec.info, "FS");
    let sor = info_f64(&rec.info, "SOR");
    let mq = info_f64(&rec.info, "MQ");
    let rp = info_f64(&rec.info, "ReadPosRankSum");
    let ic = info_f64(&rec.info, "InbreedingCoeff");
    kv(
        "rust_final",
        format!(
            "FS={fs} SOR={sor} MQ={mq} ReadPosRankSum={rp} InbreedingCoeff={ic} QUAL={:?}",
            rec.quality
        ),
    );
    assert!(keys.contains("ReadPosRankSum"));
    assert!(
        !keys.contains("InbreedingCoeff"),
        "6R.178: n=1 omits InbreedingCoeff"
    );
    assert!(fs < 0.02, "FS={fs}");
    assert!((sor - 1.179).abs() < 0.001, "SOR={sor}");
    assert!((mq - 40.25).abs() < 1e-12, "MQ={mq}");
    assert!((rec.quality.unwrap_or(0.0) - 121.84).abs() < 0.02);
}
