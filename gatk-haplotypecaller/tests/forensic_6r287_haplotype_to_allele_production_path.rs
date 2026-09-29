//! 6R.287: production-path counterfactual. Global haplotype floor, then max
//! within the live allele pools, then the existing site emitter.
//! PairHMM, genotype-likelihood construction, and rounding stay closed.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r287_haplotype_to_allele_production_path -- --nocapture --test-threads=1
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

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R287\t{key}\t{}", value.as_ref());
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

const JAVA_ALLELES: &str = include_str!("6r285_java_genotyping.tsv");
const FIRST_QNAME: &str = "HISEQ1:11:H8GV6ADXX:1:2116:18670:99941";
const JAVA_HOM_CONT: f64 = 3517.46536813194416027;
const DIAGNOSTIC_CF_CONT: f64 = 3517.46544109771957665;

struct FlagGuard;
impl Drop for FlagGuard {
    fn drop(&mut self) {
        gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r287_global_hap_floor(None, None);
    }
}

fn bits_cell(cell: &str) -> f64 {
    let hex = cell.split("bits=").nth(1).unwrap().trim();
    f64::from_bits(u64::from_str_radix(hex.trim_start_matches("0x"), 16).unwrap())
}

fn fmt_f64(x: f64) -> String {
    format!("{x:.17} bits=0x{:016x}", x.to_bits())
}

#[test]
fn forensic_6r287_haplotype_to_allele_production_path() {
    let _guard = FlagGuard;
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128);
    assert_eq!(
        HcGenotypingConfig::strict_java().stand_emit_confidence,
        DEFAULT_STAND_CALL_CONF
    );

    let gkl20 = libc_powf(10.0, -(20.0f32) / 10.0f32);
    assert_eq!(gkl20.to_bits(), Q20_GKL_POWF_F32);
    kv(
        "first_mx_primitive",
        "haplotype-to-allele max after best-4.5 floor; PairHMM closed",
    );

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
    gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r287_global_hap_floor(
        Some(TARGET),
        Some(FIRST_QNAME),
    );
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
    let emitted_pl = emitted.samples[0].pl.clone().unwrap_or_default();
    let cap = gatk_haplotypecaller::hc_genotyping_engine::take_forensic_6r287_capture()
        .expect("6R.287 capture");
    let mut java_rows = std::collections::HashMap::<(String, u16), (f64, f64, f64)>::new();
    for line in JAVA_ALLELES.lines() {
        let f: Vec<_> = line.split('\t').collect();
        if f.len() >= 7 && f[0] == "6R285" && f[1] == "allele_ll" {
            java_rows.insert(
                (f[2].to_string(), f[3].parse().unwrap()),
                (bits_cell(f[4]), bits_cell(f[5]), bits_cell(f[6])),
            );
        }
    }
    let mut matched = 0usize;
    let mut max_abs = [0.0f64; 3];
    let mut first_delta: Option<(String, u16, [f64; 3], [f64; 3])> = None;
    for row in &cap.rows {
        let Some(&(jt, jtt, jtg)) = java_rows.get(&(row.qname.clone(), row.flags)) else {
            continue;
        };
        matched += 1;
        let rust = [
            row.alleles.first().copied().unwrap_or(f64::NAN),
            row.alleles.get(1).copied().unwrap_or(f64::NAN),
            row.alleles.get(2).copied().unwrap_or(f64::NAN),
        ];
        let java = [jt, jtt, jtg];
        for i in 0..3 {
            let d = (rust[i] - java[i]).abs();
            if d > max_abs[i] {
                max_abs[i] = d;
            }
            if first_delta.is_none() && d > 0.1 {
                first_delta = Some((row.qname.clone(), row.flags, rust, java));
            }
        }
    }
    let hom_cont = cap.subset_continuous_pl.get(2).copied().unwrap_or(f64::NAN);
    let hom_pl = cap.subset_pl.get(2).copied();
    let causal = hom_pl == Some(3517) && (hom_cont - DIAGNOSTIC_CF_CONT).abs() < 1e-3;
    let classification = if causal {
        "DOWNSTREAM_HAPLOTYPE_TO_ALLELE_NORMALIZATION_CAUSAL"
    } else {
        "DOWNSTREAM_HAPLOTYPE_TO_ALLELE_PRODUCTION_PATH_NOT_JAVA_PL"
    };
    kv("n_reads", cap.n_reads.to_string());
    kv("n_haps", cap.n_haps.to_string());
    kv("alleles", cap.alleles.join(","));
    kv("ref_hap_index", cap.ref_hap_index.to_string());
    kv("n_floor_lifts", cap.n_floor_lifts.to_string());
    kv(
        "hap_allele",
        cap.hap_allele
            .iter()
            .enumerate()
            .map(|(i, a)| format!("{i}:{a}"))
            .collect::<Vec<_>>()
            .join(","),
    );
    kv(
        "allele_sums",
        cap.allele_sums
            .iter()
            .map(|v| fmt_f64(*v))
            .collect::<Vec<_>>()
            .join("\t"),
    );
    kv(
        "subset_continuous_pl",
        cap.subset_continuous_pl
            .iter()
            .map(|v| fmt_f64(*v))
            .collect::<Vec<_>>()
            .join("\t"),
    );
    kv(
        "subset_pl",
        cap.subset_pl
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(","),
    );
    kv(
        "emitted_pl",
        emitted_pl
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(","),
    );
    kv("matched_java_rows", matched.to_string());
    kv(
        "max_abs_vs_java_alleles",
        max_abs
            .iter()
            .map(|v| fmt_f64(*v))
            .collect::<Vec<_>>()
            .join("\t"),
    );
    if let Some(first) = &cap.first {
        kv("first_qname", &first.qname);
        kv("first_flags", first.flags.to_string());
        kv("first_best", fmt_f64(first.best));
        kv("first_floor", fmt_f64(first.floor));
        kv(
            "first_before",
            first
                .before
                .iter()
                .map(|v| fmt_f64(*v))
                .collect::<Vec<_>>()
                .join("\t"),
        );
        kv(
            "first_after",
            first
                .after
                .iter()
                .map(|v| fmt_f64(*v))
                .collect::<Vec<_>>()
                .join("\t"),
        );
        kv(
            "first_alleles",
            first
                .alleles
                .iter()
                .map(|v| fmt_f64(*v))
                .collect::<Vec<_>>()
                .join("\t"),
        );
    } else {
        kv("first_qname", "ABSENT");
    }
    if let Some((q, flags, rust, java)) = &first_delta {
        kv(
            "first_allele_delta",
            format!(
                "{q}\t{flags}\trust={}\tjava={}",
                rust.iter()
                    .map(|v| fmt_f64(*v))
                    .collect::<Vec<_>>()
                    .join(","),
                java.iter()
                    .map(|v| fmt_f64(*v))
                    .collect::<Vec<_>>()
                    .join(","),
            ),
        );
    }
    kv("baseline_continuous_PL", "3519.15718129777269496");
    kv("baseline_integer_PL", "3519");
    kv("counterfactual_continuous_PL", fmt_f64(hom_cont));
    kv(
        "counterfactual_integer_PL",
        hom_pl
            .map(|v| v.to_string())
            .unwrap_or_else(|| "NONE".into()),
    );
    kv("java_continuous_PL", fmt_f64(JAVA_HOM_CONT));
    kv("java_integer_PL", "3517");
    kv(
        "movement_toward_java",
        fmt_f64(3519.15718129777269496 - hom_cont),
    );
    kv("residual", fmt_f64((hom_cont - JAVA_HOM_CONT).abs()));
    kv("classification", classification);
    kv("production_change", "NONE");
    assert_eq!(cap.n_reads, 123);
    assert_eq!(cap.n_floor_lifts, 0);
    assert_eq!(cap.subset_pl, vec![570, 0, 3518]);
    assert_eq!(emitted_pl.as_slice(), &[570, 0, 3518]);
    assert_eq!(
        classification,
        "DOWNSTREAM_HAPLOTYPE_TO_ALLELE_PRODUCTION_PATH_NOT_JAVA_PL"
    );
    if cap.n_haps == 24 && cap.alleles == ["T", "TTTG", "TGTTTG"] {
        for i in 0..14 {
            assert_eq!(cap.hap_allele[i], "T", "hap {i}");
        }
        for i in 14..19 {
            assert_eq!(cap.hap_allele[i], "TTTG", "hap {i}");
        }
        for i in 19..24 {
            assert_eq!(cap.hap_allele[i], "TGTTTG", "hap {i}");
        }
        assert_eq!(cap.ref_hap_index, 4);
    }
}
