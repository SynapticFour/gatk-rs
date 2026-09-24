//! 6R.257 live: under frozen Java 174-mers + Java clip, PairHMM primitive
//! input arrays match Java 4.4 `modifyReadQualities`+GCP on the 123-read
//! `20:29455649 T/TGTTTG` object. IQ/DQ Q6 floor never fires. Remaining
//! difference is kernel configuration (GKL vs NEON f64), not trim/clip.
//! Proof-only. Skipped unless `HOLDOUT_6R257=1`.
//!
//! ```text
//! HOLDOUT_6R257=1 cargo test -p gatk-haplotypecaller --test holdout_6r257_pairhmm_input_plane -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::assembly_region_finalize::{
    clip_finalized_reads_to_region, finalize_region_reads_for_assembly,
    gatk_min_tail_quality_for_assembly,
};
use gatk_haplotypecaller::hc_genotyping_engine::{
    take_colocated_merge_numerics, HcGenotypingConfig,
};
use gatk_haplotypecaller::pairhmm_qual::MIN_USABLE_Q_SCORE;
use gatk_haplotypecaller::pcr_error_model::{apply_pcr_error_model, PcrErrorModel};
use gatk_haplotypecaller::read_threading_assembler::DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;
use gatk_haplotypecaller::{
    begin_hap_list_observe, call_disposition, flatten_assembly_regions,
    indel_gop_from_optional_tag, take_hap_list_trim_span, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, HcLikelihoodEngineConfig, PairHmmBackend, ReadFilterParams,
    WalkerTraversalConfig, DEFAULT_STAND_CALL_CONF, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use rust_htslib::bam::record::Aux;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/java.vcf";
const TARGET: u64 = 29_455_649;
const TARGET_REF: &str = "T";
const TARGET_ALT: &str = "TGTTTG";
const CLOSED_GC: u64 = 29_455_314;
const FROZEN_IDX: [usize; 5] = [19, 20, 21, 22, 23];
const RUST_FNV: [&str; 5] = [
    "55012fcf3b430591",
    "341b2e070ccb5846",
    "6a65e4c02733c2ed",
    "fa07750bf228b3c2",
    "79451c576721a729",
];
const JAVA_FNV: [&str; 5] = [
    "fcd72c6ce600dd16",
    "c1a4e8204522f645",
    "eb03271fa7548f26",
    "7dc2a8ae5da116e1",
    "343e7c443c1d9732",
];
const RUST_CLIP: (u64, u64) = (29_455_569, 29_455_724);
const JAVA_CLIP: (u64, u64) = (29_455_560, 29_455_728);
const JAVA_TSV: &str = include_str!("forensic_6r252_java_hap_trim.tsv");

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R257\t{key}\t{}", value.as_ref());
}

fn fnv1a64_hex(bases: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bases {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{h:016x}")
}

fn java_seq(hash: &str) -> Vec<u8> {
    for line in JAVA_TSV.lines() {
        if line.contains("stage=trimmed") && line.contains(&format!("hash={hash}")) {
            if let Some(seq) = line.split("\tseq=").nth(1) {
                return seq.as_bytes().to_vec();
            }
        }
    }
    panic!("missing Java trimmed seq for {hash}");
}

fn bam_indel_phred(rec: &rust_htslib::bam::Record, tag: &[u8]) -> Option<Vec<u8>> {
    match rec.aux(tag) {
        Ok(Aux::String(s)) => Some(s.bytes().map(|b| b.saturating_sub(33)).collect()),
        _ => None,
    }
}

fn clip_map(
    finalized: &[rust_htslib::bam::Record],
    region: &gatk_haplotypecaller::assembly_region_iterator::AssemblyRegion,
    start: u64,
    end: u64,
) -> HashMap<(Vec<u8>, u16), rust_htslib::bam::Record> {
    let mut clip_region = region.clone();
    clip_region.extended_start = GenomePosition::new_1based(start);
    clip_region.extended_end = GenomePosition::new_1based(end);
    let mut clipped = clip_finalized_reads_to_region(finalized, &clip_region);
    clipped.retain(|r| r.seq().len() >= 10);
    let mut by_id = HashMap::new();
    for rec in clipped {
        by_id.insert((rec.qname().to_vec(), rec.flags()), rec);
    }
    by_id
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
fn holdout_6r257_pairhmm_input_plane() {
    if std::env::var("HOLDOUT_6R257").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R257=1");
        return;
    }
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128);
    assert_eq!(
        HcGenotypingConfig::strict_java().stand_emit_confidence,
        DEFAULT_STAND_CALL_CONF
    );
    let cfg = HcLikelihoodEngineConfig::gatk_haplotype_caller_production();
    assert_eq!(cfg.resolved_pair_hmm_backend(), PairHmmBackend::NeonF64);
    assert_eq!(cfg.pcr_error_model, PcrErrorModel::Conservative);

    let root = repo_root();
    assert!(vcf_has(
        &root.join(JAVA_VCF_REL),
        TARGET,
        TARGET_REF,
        TARGET_ALT
    ));
    assert!(!vcf_has(&root.join(JAVA_VCF_REL), CLOSED_GC, "G", "C"));

    let java_haps: Vec<Vec<u8>> = JAVA_FNV.iter().map(|h| java_seq(h)).collect();
    for (i, seq) in java_haps.iter().enumerate() {
        assert_eq!(fnv1a64_hex(seq), JAVA_FNV[i]);
        assert_eq!(seq.len(), 174);
    }

    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, INTERVAL).expect("interval");
    begin_hap_list_observe();
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
    let trim = take_hap_list_trim_span().expect("trim");
    assert_eq!((trim.trim_start, trim.trim_end), RUST_CLIP);
    let snap = take_colocated_merge_numerics()
        .into_iter()
        .find(|s| s.loc == TARGET)
        .expect("snap");
    assert_eq!(snap.n_reads, 123);
    for (i, &idx) in FROZEN_IDX.iter().enumerate() {
        let rust = &outcome.assembly.haplotypes[idx].bases;
        assert_eq!(fnv1a64_hex(rust), RUST_FNV[i]);
        assert_eq!(&java_haps[i][9..java_haps[i].len() - 4], rust.as_slice());
    }
    let recs =
        try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let emitted = recs
        .iter()
        .find(|r| {
            r.position == TARGET
                && r.reference == TARGET_REF
                && r.alternate.iter().any(|a| a == TARGET_ALT)
        })
        .expect("emit");
    assert_eq!(
        emitted.samples[0].pl.as_deref(),
        Some(&[570, 0, 3518][..]),
        "production PL must remain 3518"
    );

    let finalized = finalize_region_reads_for_assembly(
        &region.reads,
        region,
        true,
        gatk_min_tail_quality_for_assembly(10),
        false,
    );
    let java_map = clip_map(&finalized, region, JAVA_CLIP.0, JAVA_CLIP.1);
    let mut n_iq_lt6 = 0usize;
    let mut n_dq_lt6 = 0usize;
    let mut n_iq_floor_noop = 0usize;
    for (ri, _) in snap.ad_row_read_index.iter().enumerate() {
        let key = (
            snap.ad_row_qname[ri].as_bytes().to_vec(),
            snap.ad_row_flags[ri],
        );
        let rec = java_map.get(&key).expect("Java-clip read");
        let n = rec.seq().len();
        let mut iq =
            indel_gop_from_optional_tag(bam_indel_phred(rec, b"BI").as_deref(), n).unwrap();
        let mut dq =
            indel_gop_from_optional_tag(bam_indel_phred(rec, b"BD").as_deref(), n).unwrap();
        apply_pcr_error_model(&rec.seq().as_bytes(), &mut iq, &mut dq, cfg.pcr_error_model);
        n_iq_lt6 += iq.iter().filter(|&&q| q < MIN_USABLE_Q_SCORE).count();
        n_dq_lt6 += dq.iter().filter(|&&q| q < MIN_USABLE_Q_SCORE).count();
        if iq.iter().all(|&q| q >= MIN_USABLE_Q_SCORE)
            && dq.iter().all(|&q| q >= MIN_USABLE_Q_SCORE)
        {
            n_iq_floor_noop += 1;
        }
    }
    assert_eq!(n_iq_lt6, 0);
    assert_eq!(n_dq_lt6, 0);
    assert_eq!(n_iq_floor_noop, 123);
    kv("production_integer_PL", "3518");
    kv("n_iq_bytes_lt_q6_after_pcr", "0");
    kv("n_dq_bytes_lt_q6_after_pcr", "0");
    kv("classification", "PAIRHMM_CONFIGURATION_DIVERGENCE");
    kv("production_change", "NONE");
}
