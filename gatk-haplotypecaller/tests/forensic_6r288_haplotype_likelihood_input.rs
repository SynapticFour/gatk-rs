//! 6R.288: trace one read's haplotype likelihoods from the PairHMM return
//! to `try_genotype_colocated_snp_indel_merge`, then substitute Java's
//! 24-column row at that boundary.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r288_haplotype_likelihood_input -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::hc_genotyping_engine::HcGenotypingConfig;
use gatk_haplotypecaller::read_threading_assembler::DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;
use gatk_haplotypecaller::{
    begin_hap_list_observe, begin_likelihood_pipeline_observe, call_disposition,
    flatten_assembly_regions, take_hap_list_trim_span, take_likelihood_pipeline_cells,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, HaplotypeCallerEngine, HcLikelihoodEngineConfig, LikelihoodPipelineCell,
    PairHmmBackend, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_CALL_CONF,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::collections::BTreeMap;
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
const JAVA_HAP_ROW: &str = include_str!("6r286_haplotype_allele_map.tsv");
const Q20_GKL_POWF_F32: u32 = 0x3c23d70a;
const FIRST_QNAME: &str = "HISEQ1:11:H8GV6ADXX:1:2116:18670:99941";
const FIRST_FLAGS: u16 = 99;
const SHOW: [usize; 9] = [0, 1, 4, 9, 19, 20, 21, 22, 23];
const RUST_H0_BITS: u64 = 0xc004aa3284709780;
const JAVA_H0_BITS: u64 = 0xc00d9faa00000000;
const JAVA_HOM_CONT: f64 = 3517.46536813194416027;
const EMITTER_BASELINE_CONT: f64 = 3517.51626740832807627;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R288\t{key}\t{}", value.as_ref());
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
        let mut c = line.split('\t');
        let Some(chrom) = c.next() else { continue };
        let Some(p) = c.next() else { continue };
        let Some(_id) = c.next() else { continue };
        let Some(reference) = c.next() else { continue };
        let Some(alts) = c.next() else { continue };
        if chrom == "20"
            && p == pos.to_string()
            && reference == r
            && alts.split(',').any(|x| x == a)
        {
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

fn fmt_f64(x: f64) -> String {
    format!("{x:.17} bits=0x{:016x}", x.to_bits())
}

fn java_hap_row() -> Vec<f64> {
    for line in JAVA_HAP_ROW.lines() {
        let f: Vec<_> = line.split('\t').collect();
        if f.len() >= 28
            && f[0] == "6R286"
            && f[1] == "hap_row"
            && f[2] == FIRST_QNAME
            && f[3] == "99"
        {
            return f[4..28]
                .iter()
                .map(|hex| f64::from_bits(u64::from_str_radix(hex, 16).unwrap()))
                .collect();
        }
    }
    panic!("missing Java hap_row for {FIRST_QNAME}");
}

struct StageSnap {
    seq: u32,
    stage: &'static str,
    n_reads: usize,
    n_haps: usize,
    read_index: usize,
    cells: BTreeMap<usize, (u64, f64)>,
}

fn stages_for(cells: &[LikelihoodPipelineCell]) -> Vec<StageSnap> {
    let mut order: Vec<(u32, &'static str)> = Vec::new();
    let mut grouped: BTreeMap<(u32, &'static str), StageSnap> = BTreeMap::new();
    for cell in cells {
        if cell.qname != FIRST_QNAME || cell.flags != FIRST_FLAGS {
            continue;
        }
        let key = (cell.seq, cell.stage);
        if !grouped.contains_key(&key) {
            order.push(key);
            grouped.insert(
                key,
                StageSnap {
                    seq: cell.seq,
                    stage: cell.stage,
                    n_reads: cell.n_reads,
                    n_haps: cell.n_haps,
                    read_index: cell.read_index,
                    cells: BTreeMap::new(),
                },
            );
        }
        grouped
            .get_mut(&key)
            .unwrap()
            .cells
            .insert(cell.hap_index, (cell.hap_fnv, cell.log10_likelihood));
    }
    order
        .into_iter()
        .filter_map(|k| grouped.remove(&k))
        .collect()
}

fn best_of(stage: &StageSnap) -> (usize, f64) {
    stage
        .cells
        .iter()
        .filter(|(_, (_, v))| v.is_finite())
        .max_by(|a, b| a.1 .1.total_cmp(&b.1 .1))
        .map(|(&i, &(_, v))| (i, v))
        .unwrap_or((usize::MAX, f64::NAN))
}

fn cell_bits(stage: &StageSnap, hap: usize) -> Option<u64> {
    stage.cells.get(&hap).map(|(_, v)| v.to_bits())
}

struct FlagGuard;
impl Drop for FlagGuard {
    fn drop(&mut self) {
        gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r287_global_hap_floor(None, None);
        gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r288_haplotype_substitute(
            None, None, 0, None,
        );
    }
}

#[test]
fn forensic_6r288_haplotype_likelihood_input() {
    let _guard = FlagGuard;
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128);
    assert_eq!(
        HcGenotypingConfig::strict_java().stand_emit_confidence,
        DEFAULT_STAND_CALL_CONF
    );
    let gkl20 = libc_powf(10.0, -(20.0f32) / 10.0f32);
    assert_eq!(gkl20.to_bits(), Q20_GKL_POWF_F32);

    let java_row = java_hap_row();
    assert_eq!(java_row.len(), 24);
    assert_eq!(java_row[0].to_bits(), JAVA_H0_BITS);
    assert_eq!(java_row[1].to_bits(), JAVA_H0_BITS);

    let root = repo_root();
    let java_vcf = root.join(JAVA_VCF_REL);
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !java_vcf.is_file() || !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }
    assert!(vcf_has(&java_vcf, TARGET, TARGET_REF, TARGET_ALT));
    assert!(!vcf_has(&java_vcf, CLOSED_GC, "G", "C"));
    let java_haps: Vec<Vec<u8>> = JAVA_FNV.iter().map(|h| java_seq(h)).collect();
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
    gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r287_global_hap_floor(None, None);
    gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r288_haplotype_substitute(
        Some(TARGET),
        Some(FIRST_QNAME),
        FIRST_FLAGS,
        Some(java_row.clone()),
    );
    begin_likelihood_pipeline_observe();
    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let cells = take_likelihood_pipeline_cells();
    let snaps = stages_for(&cells);
    let trim = take_hap_list_trim_span().expect("trim");
    assert_eq!((trim.trim_start, trim.trim_end), RUST_CLIP);
    for (i, &idx) in FROZEN_IDX.iter().enumerate() {
        let rust = &outcome.assembly.haplotypes[idx].bases;
        assert_eq!(fnv1a64_hex(rust), RUST_FNV[i]);
        assert_eq!(&java_haps[i][9..java_haps[i].len() - 4], rust.as_slice());
    }
    assert_eq!(
        HcLikelihoodEngineConfig::gatk_haplotype_caller_production().resolved_pair_hmm_backend(),
        PairHmmBackend::NeonF64
    );
    let report = gatk_haplotypecaller::hc_genotyping_engine::take_forensic_6r288_report()
        .expect("6R.288 report");
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

    kv("qname", FIRST_QNAME);
    kv("flags", FIRST_FLAGS.to_string());
    kv("merge_read_index", report.read_index.to_string());
    kv("n_stages", snaps.len().to_string());
    let frozen_fnv = u64::from_str_radix(RUST_FNV[0], 16).unwrap();
    for stage in &snaps {
        let (best_i, best) = best_of(stage);
        let cols = SHOW
            .iter()
            .map(|&i| match stage.cells.get(&i) {
                Some((_, v)) => format!("H{i}={}", fmt_f64(*v)),
                None => format!("H{i}=ABSENT"),
            })
            .collect::<Vec<_>>()
            .join("\t");
        let same_final = stage
            .cells
            .get(&19)
            .is_some_and(|(fnv, _)| *fnv == frozen_fnv);
        kv(
            "checkpoint",
            format!(
                "seq={}\tstage={}\tn_reads={}\tn_haps={}\tread_index={}\tbest_hap={best_i}\tbest={}\tsame_final_haps={same_final}\t{cols}",
                stage.seq,
                stage.stage,
                stage.n_reads,
                stage.n_haps,
                stage.read_index,
                fmt_f64(best),
            ),
        );
    }
    let (best_i, best) = {
        let mut best_i = 0usize;
        let mut best = f64::NEG_INFINITY;
        for (i, v) in report.before.iter().enumerate() {
            if v.is_finite() && *v > best {
                best = *v;
                best_i = i;
            }
        }
        (best_i, best)
    };
    kv(
        "checkpoint",
        format!(
            "seq=merge\tstage=merge_entry\tn_reads=123\tn_haps={}\tread_index={}\tbest_hap={best_i}\tbest={}\tsame_final_haps=true\t{}",
            report.before.len(),
            report.read_index,
            fmt_f64(best),
            SHOW.iter()
                .map(|&i| report
                    .before
                    .get(i)
                    .map(|v| format!("H{i}={}", fmt_f64(*v)))
                    .unwrap_or_else(|| format!("H{i}=ABSENT")))
                .collect::<Vec<_>>()
                .join("\t"),
        ),
    );

    let final_stages: Vec<&StageSnap> = snaps
        .iter()
        .filter(|s| s.cells.get(&19).is_some_and(|(fnv, _)| *fnv == frozen_fnv))
        .collect();
    let scoring = final_stages
        .iter()
        .rev()
        .find(|s| s.stage == "post_kernel" || s.stage == "refresh");
    let normalize = final_stages.iter().rev().find(|s| s.stage == "normalize");
    let filter = final_stages.iter().rev().find(|s| s.stage == "filter");
    let score_h0 = scoring.and_then(|s| cell_bits(s, 0));
    let norm_h0 = normalize.and_then(|s| cell_bits(s, 0));
    let filter_h0 = filter.and_then(|s| cell_bits(s, 0));
    let merge_h0 = report.before.first().map(|v| v.to_bits());
    kv(
        "h0_bits",
        format!(
            "score={score_h0:?}\tnormalize={norm_h0:?}\tfilter={filter_h0:?}\tmerge={merge_h0:?}\tjava={JAVA_H0_BITS:#x}"
        ),
    );
    let chain_equal = score_h0 == norm_h0 && norm_h0 == filter_h0 && filter_h0 == merge_h0;
    let h0_is_rust = merge_h0 == Some(RUST_H0_BITS);
    let t_pool_max = report.before[..report.before.len().min(14)]
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .fold(f64::NEG_INFINITY, f64::max);
    kv("t_pool_max", fmt_f64(t_pool_max));
    kv(
        "h4",
        report
            .before
            .get(4)
            .map(|v| fmt_f64(*v))
            .unwrap_or_else(|| "ABSENT".into()),
    );
    kv("chain_equal", chain_equal.to_string());
    let hom = report
        .subset_continuous_pl
        .get(2)
        .copied()
        .unwrap_or(f64::NAN);
    let gap_before = (EMITTER_BASELINE_CONT - JAVA_HOM_CONT).abs();
    let gap_after = (hom - JAVA_HOM_CONT).abs();
    let toward = gap_before - gap_after;
    let classification = if chain_equal && h0_is_rust && toward < 1e-3 {
        "HAPLOTYPE_LIKELIHOOD_PAIRHMM_RETURN_NOT_CAUSAL"
    } else if chain_equal && h0_is_rust && toward >= 1e-3 {
        "HAPLOTYPE_LIKELIHOOD_INPUT_CAUSAL"
    } else {
        "HAPLOTYPE_LIKELIHOOD_CHECKPOINT_NOT_ISOLATED"
    };
    kv("classification", classification);
    kv("counterfactual_continuous_PL", fmt_f64(hom));
    kv(
        "counterfactual_integer_PL",
        report
            .subset_pl
            .get(2)
            .map(|v| v.to_string())
            .unwrap_or_else(|| "NONE".into()),
    );
    kv(
        "emitted_pl",
        emitted_pl
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(","),
    );
    kv("movement_toward_java", fmt_f64(toward));
    kv("residual", fmt_f64(gap_after));
    kv(
        "delta_h0",
        fmt_f64(f64::from_bits(RUST_H0_BITS) - f64::from_bits(JAVA_H0_BITS)),
    );
    assert!(chain_equal, "H0 changed between PairHMM return and merge");
    assert!(h0_is_rust, "merge H0 is not the known production value");
    assert_eq!(report.before[0].to_bits(), report.before[1].to_bits());
    assert_eq!(t_pool_max.to_bits(), report.before[0].to_bits());
    assert_ne!(report.before[4].to_bits(), report.before[0].to_bits());
    assert_eq!(
        classification,
        "HAPLOTYPE_LIKELIHOOD_PAIRHMM_RETURN_NOT_CAUSAL"
    );
}
