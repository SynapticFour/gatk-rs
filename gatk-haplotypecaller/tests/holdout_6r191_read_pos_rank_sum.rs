//! 6R.191 live: `2:92305716 A/C` omits ReadPosRankSum from annotation likelihoods.
//! FORMAT stays GT=1/1 AD=0,3 DP=3 GQ=9 PL=130,9,0 QUAL=116.84. Annotation n=4.
//! Skipped unless `HOLDOUT_6R191=1`.
//!
//! ```text
//! HOLDOUT_6R191=1 cargo test -p gatk-haplotypecaller --test holdout_6r191_read_pos_rank_sum -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::variant_site_hc_annotations::{
    read_pos_rank_sum_quals_from_likelihoods, HcStrandBiasLikelihoods,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, HcGenotypingConfig, ReadFilterParams, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const GAP_INTERVAL: &str = "2:92305500-92305850";
const TG_INTERVAL: &str = "2:92307200-92307550";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_305_716;
const CLOSED_SNP: u64 = 92_305_634;
const CLOSED_AG: u64 = 92_305_635;
const CLOSED_TG: u64 = 92_307_333;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R191\t{key}\t{}", value.as_ref());
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

fn info_has(info: &[InfoValue], key: &str) -> bool {
    info.iter().any(|v| match v {
        InfoValue::Flag(k)
        | InfoValue::Integer(k, _)
        | InfoValue::Float(k, _)
        | InfoValue::String(k, _)
        | InfoValue::Character(k, _) => k == key,
    })
}

#[test]
fn holdout_6r191_read_pos_rank_sum() {
    if std::env::var("HOLDOUT_6R191").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R191=1");
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
        "ReadPosRankSum from annotation_likelihoods fillQualsFromLikelihood",
    );
    kv("target", "2:92305716 A/C");
    kv("classification", "A — WRONG SOURCE OBJECT (closed)");

    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, GAP_INTERVAL).expect("interval");
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
    let covering_closed = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_AG
                && r.end.get() >= CLOSED_SNP
        })
        .expect("covering closed");
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
    let closed_outcome = HaplotypeCallerEngine::call_region(
        covering_closed,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("closed call")
    .expect("closed outcome");
    let outcome = HaplotypeCallerEngine::call_region(
        covering,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("outcome");

    let closed_ag = closed_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_AG)
                && c.event.ref_allele == "A"
                && c.event.alt_allele == "G"
        })
        .expect("A/G");
    let target = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(TARGET)
                && c.event.ref_allele == "A"
                && c.event.alt_allele == "C"
        })
        .expect("A/C");
    let closed_ag_ann: BTreeSet<usize> = closed_ag
        .annotation_likelihoods
        .iter()
        .map(|c| c.read_index.get())
        .collect();
    let target_ann: BTreeSet<usize> = target
        .annotation_likelihoods
        .iter()
        .map(|c| c.read_index.get())
        .collect();
    assert_eq!(closed_ag_ann.len(), 3);
    assert_eq!(target_ann.len(), 4);
    assert_eq!(target.genotype.format.ad_as_i32(), vec![0, 3]);
    assert_eq!(target.genotype.format.dp.as_i32(), 3);
    assert_eq!(target.genotype.format.pl_as_i32(), vec![130, 9, 0]);

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
        likelihoods: &target.annotation_likelihoods,
        haplotypes: &outcome.assembly.haplotypes,
        contig: &covering.contig,
        ref_bytes: outcome.assembly.reference_bases(),
        pad_start_1based: apply_pad,
        full_ref_bytes: full_ref,
        full_pad_1based: full_pad,
        max_mnp_distance: outcome.assembly.max_mnp_distance(),
        emit_spanning_dels: !cfg.disable_spanning_event_genotyping,
    };
    let (refs, alts) = read_pos_rank_sum_quals_from_likelihoods(&evidence, TARGET, "A", "C");
    assert_eq!(refs.len(), 0);
    assert_eq!(alts.len(), 3);

    let closed_emitted = try_emit_call_region_variants(
        covering_closed,
        &closed_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("closed emit");
    let emitted =
        try_emit_call_region_variants(covering, &outcome, "NA12878", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let closed_ag_rec = closed_emitted
        .iter()
        .find(|r| r.position == CLOSED_AG && r.reference == "A")
        .expect("A/G emit");
    let target_rec = emitted
        .iter()
        .find(|r| r.position == TARGET && r.reference == "A")
        .expect("A/C emit");
    let sample = target_rec.samples.first().expect("sample");
    assert_eq!(
        sample.gt.as_ref().map(|g| g.to_string()).as_deref(),
        Some("1/1")
    );
    assert_eq!(sample.ad.as_deref(), Some(&[0u32, 3][..]));
    assert_eq!(sample.dp.map(|d| d as i32), Some(3));
    assert!((target_rec.quality.unwrap_or(0.0) - 116.84).abs() < 0.05);
    assert_eq!(info_i32(&closed_ag_rec.info, "DP"), Some(3));
    let ag_mq = info_f64(&closed_ag_rec.info, "MQ").unwrap_or(-1.0);
    let ag_sor = info_f64(&closed_ag_rec.info, "SOR").unwrap_or(-1.0);
    assert!((ag_mq - 41.96).abs() < 0.005);
    assert!((ag_sor - 0.693).abs() < 0.002);
    assert!(!info_has(&target_rec.info, "ReadPosRankSum"));
    kv("format", "GT=1/1 AD=0,3 DP=3 GQ=9 PL=130,9,0 QUAL=116.84");
    kv("info", "DP=4 MQ=35.15 ReadPosRankSum=absent");

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
