//! 6R.244 live: colocated merge at `20:29455649 T/TGTTTG`. 6R.246 closed
//! the empty-annotation gap: unique evidence 123, INFO DP=123.
//! Skipped unless `HOLDOUT_6R244=1`.
//!
//! ```text
//! HOLDOUT_6R244=1 cargo test -p gatk-haplotypecaller --test holdout_6r244_colocated_merge_annotation -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::hc_genotyping_engine::{
    take_colocated_merge_numerics, HcGenotypingConfig,
};
use gatk_haplotypecaller::read_threading_assembler::DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_CALL_CONF,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/java.vcf";
const TARGET: u64 = 29_455_649;
const TARGET_REF: &str = "T";
const TARGET_ALT: &str = "TGTTTG";
const CLOSED_GC: u64 = 29_455_314;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R244\t{key}\t{}", value.as_ref());
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

fn vcf_has(path: &Path, pos: u64, r: &str, a: &str) -> bool {
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

#[test]
fn holdout_6r244_colocated_merge_annotation() {
    if std::env::var("HOLDOUT_6R244").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R244=1");
        return;
    }
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128);
    assert_eq!(
        HcGenotypingConfig::strict_java().stand_emit_confidence,
        DEFAULT_STAND_CALL_CONF
    );
    let root = repo_root();
    let java_vcf = root.join(JAVA_VCF_REL);
    assert!(vcf_has(&java_vcf, TARGET, TARGET_REF, TARGET_ALT));
    assert!(!vcf_has(&java_vcf, CLOSED_GC, "G", "C"));

    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
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
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("ActiveFull");
    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let at_loc: Vec<_> = outcome
        .genotyped_calls
        .iter()
        .filter(|c| c.event.start_1based.get() == TARGET)
        .collect();
    assert_eq!(at_loc.len(), 1, "SiteScore did not add a second call");
    let call = at_loc[0];
    assert_eq!(call.event.ref_allele, TARGET_REF);
    assert_eq!(call.event.alt_allele, TARGET_ALT);
    assert!(call.post_merge_unused_alt_subset);
    let unique: BTreeSet<usize> = call
        .annotation_likelihoods
        .iter()
        .map(|c| c.read_index.get())
        .collect();
    assert_eq!(unique.len(), 123, "6R.246 attached the merge subset");
    assert_eq!(call.annotation_likelihoods.len(), 2952);
    assert_eq!(call.genotype.format.ad_as_i32(), vec![88, 22]);
    let snap = take_colocated_merge_numerics()
        .into_iter()
        .find(|s| s.loc == TARGET)
        .expect("merge snapshot");
    assert_eq!(snap.n_overlap_before_qname_dedupe, 123);
    assert_eq!(snap.n_pairhmm_reads, 230);
    let recs =
        try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let emitted = recs
        .iter()
        .find(|r| r.position == TARGET && r.reference == TARGET_REF)
        .expect("emit");
    assert_eq!(info_i32(&emitted.info, "DP"), Some(123));
    kv(
        "classification",
        "ANNOTATION_LIKELIHOOD_LIFECYCLE_DIVERGENCE",
    );
    kv(
        "first_divergent_operation",
        "try_genotype_colocated_snp_indel_merge annotation_likelihoods=subset.into_owned() (6R.246)",
    );
    kv("production_change", "6R.246 subset.into_owned()");
    kv("next_arrow", "closed by 6R.246");
}
