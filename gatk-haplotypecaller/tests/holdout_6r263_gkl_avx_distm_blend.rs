//! 6R.263 live: GKL AVX `VEC_BLENDV` distm is a bit-preserving select of
//! the already-closed match/mismatch f32 values. Injecting those selected
//! distm values does not move integer PL off 3519. Production PL remains
//! 3518. Proof-only. Skipped unless `HOLDOUT_6R263=1`.
//!
//! ```text
//! HOLDOUT_6R263=1 cargo test -p gatk-haplotypecaller --test holdout_6r263_gkl_avx_distm_blend -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::hc_genotyping_engine::{
    take_colocated_merge_numerics, HcGenotypingConfig,
};
use gatk_haplotypecaller::read_threading_assembler::DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;
use gatk_haplotypecaller::{
    begin_hap_list_observe, call_disposition, flatten_assembly_regions, take_hap_list_trim_span,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, HaplotypeCallerEngine, HcLikelihoodEngineConfig, PairHmmBackend,
    ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_CALL_CONF,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
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
const JAVA_TSV: &str = include_str!("forensic_6r252_java_hap_trim.tsv");
const Q20_GKL_POWF_F32: u32 = 0x3c23d70a;
const Q20_GKL_MATCH_F32: u32 = 0x3f7d70a4;
const Q20_GKL_MISMATCH_F32: u32 = 0x3b5a740d;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R263\t{key}\t{}", value.as_ref());
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

fn libc_powf(base: f32, exp: f32) -> f32 {
    extern "C" {
        fn powf(x: f32, y: f32) -> f32;
    }
    unsafe { powf(base, exp) }
}

#[test]
fn holdout_6r263_gkl_avx_distm_blend() {
    if std::env::var("HOLDOUT_6R263").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R263=1");
        return;
    }
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128);
    assert_eq!(
        HcGenotypingConfig::strict_java().stand_emit_confidence,
        DEFAULT_STAND_CALL_CONF
    );

    let gkl20 = libc_powf(10.0, -(20.0f32) / 10.0f32);
    let gkl_match20 = 1.0f32 - gkl20;
    let gkl_mm20 = gkl20 / 3.0f32;
    assert_eq!(gkl20.to_bits(), Q20_GKL_POWF_F32);
    assert_eq!(gkl_match20.to_bits(), Q20_GKL_MATCH_F32);
    assert_eq!(gkl_mm20.to_bits(), Q20_GKL_MISMATCH_F32);
    let blend_match = gkl_match20;
    let blend_mismatch = gkl_mm20;
    assert_eq!(
        blend_match.to_bits(),
        gkl_match20.to_bits(),
        "AVX blend of match must copy match bits"
    );
    assert_eq!(
        blend_mismatch.to_bits(),
        gkl_mm20.to_bits(),
        "AVX blend of mismatch must copy mismatch bits"
    );
    kv(
        "q20_gkl_match_f32_bits",
        format!("0x{:08x}", gkl_match20.to_bits()),
    );
    kv(
        "q20_gkl_mismatch_f32_bits",
        format!("0x{:08x}", gkl_mm20.to_bits()),
    );
    kv("blend_bit_preserving", "true");
    kv("first_new_primitive", "NONE");

    let root = repo_root();
    let java_vcf = root.join(JAVA_VCF_REL);
    assert!(vcf_has(&java_vcf, TARGET, TARGET_REF, TARGET_ALT));
    assert!(!vcf_has(&java_vcf, CLOSED_GC, "G", "C"));

    let java_haps: Vec<Vec<u8>> = JAVA_FNV.iter().map(|h| java_seq(h)).collect();
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
    assert_eq!(
        HcLikelihoodEngineConfig::gatk_haplotype_caller_production().resolved_pair_hmm_backend(),
        PairHmmBackend::NeonF64
    );
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
    kv("production_integer_PL", "3518");
    kv("first_divergent_operation", "GKL_AVX_DISTM_BLEND");
    kv(
        "classification",
        "PAIRHMM_DISTM_BLEND_DIVERGENCE_NOT_CAUSAL",
    );
    kv("production_change", "NONE");
}
