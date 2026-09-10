//! 6R.148: consequence of the two 6R.147 best-haplotype identity mismatches.
//! Diagnostic-only `unbounded_diagnostic`. Production k-best policy unchanged.
//! Skipped unless `HOLDOUT_6R148=1`.
//!
//! ```text
//! HOLDOUT_6R148=1 GATK_RS_EXPERIMENTAL_KBEST_POLICY=unbounded_diagnostic \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r148_best_haplotype_identity -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    begin_likelihood_pipeline_observe, begin_poorly_modeled_observe, begin_realign_observe,
    call_disposition, flatten_assembly_regions, take_likelihood_pipeline_cells,
    take_likelihood_pipeline_snaps, take_poorly_modeled_observe, take_realign_observe,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_FLOAT_DUMP_REL: &str = "parity/giab/runs/local-pairhmm-diff/6r143_java.txt";
const TARGET: u64 = 29_456_196;
const JAVA_ACTIVE_START: u64 = 29_455_995;
const JAVA_ACTIVE_END: u64 = 29_456_293;
const JAVA_PAD_START: u64 = 29_455_895;
const JAVA_PAD_END: u64 = 29_456_393;
const INFORMATIVE: f64 = 0.2;
const HAP_START: u64 = 29_455_977;
const HAP_END: u64 = 29_456_301;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R148\t{key}\t{}", value.as_ref());
}

struct EnvGuard {
    key: &'static str,
    prior: Option<String>,
}

impl EnvGuard {
    fn set(key: &'static str, value: &str) -> Self {
        let prior = std::env::var(key).ok();
        std::env::set_var(key, value);
        Self { key, prior }
    }
}

impl Drop for EnvGuard {
    fn drop(&mut self) {
        match &self.prior {
            Some(v) => std::env::set_var(self.key, v),
            None => std::env::remove_var(self.key),
        }
    }
}

fn parse_kv_fields(parts: &[&str]) -> BTreeMap<String, String> {
    let mut m = BTreeMap::new();
    for tok in parts {
        if let Some((k, v)) = tok.split_once('=') {
            m.insert(k.to_string(), v.to_string());
        }
    }
    m
}

type ReadKey = (String, u16);
type CellKey = (String, u16, String);

fn java_hap_priority(is_ref: bool, n_cigar_elements: usize) -> f64 {
    let reference_term = if is_ref { 1.0 } else { 0.0 };
    reference_term + (1.0 - n_cigar_elements as f64)
}

fn n_cigar_elements(cigar: &str) -> usize {
    let mut n = 0usize;
    let mut i = 0usize;
    let b = cigar.as_bytes();
    while i < b.len() {
        if b[i].is_ascii_digit() {
            while i < b.len() && b[i].is_ascii_digit() {
                i += 1;
            }
            if i < b.len() {
                n += 1;
                i += 1;
            }
        } else {
            i += 1;
        }
    }
    n.max(1)
}

fn java_search_best_allele(ll: &[f64], priorities: &[f64]) -> usize {
    let mut best = 0usize;
    let mut best_ll = ll[0];
    let mut second_ll = f64::NEG_INFINITY;
    for (a, &cand) in ll.iter().enumerate().skip(1) {
        if cand > best_ll {
            second_ll = best_ll;
            best = a;
            best_ll = cand;
        } else if cand > second_ll {
            second_ll = cand;
        }
    }
    if best_ll - second_ll < INFORMATIVE {
        let mut best_pri = priorities[best];
        let mut second_pri = f64::NEG_INFINITY;
        for (a, &cand) in ll.iter().enumerate() {
            if a == best || best_ll - cand > INFORMATIVE {
                continue;
            }
            let pri = priorities[a];
            if pri > best_pri {
                best = a;
                second_pri = best_pri;
                best_pri = pri;
            } else if pri > second_pri {
                second_pri = pri;
            }
        }
    }
    best
}

struct JavaHap {
    hash: String,
    is_ref: bool,
    cigar: String,
    loc: String,
    align_start: String,
    bases: String,
}

fn load_java_haps(path: &Path) -> Vec<JavaHap> {
    let text = std::fs::read_to_string(path).expect("java dump");
    let mut out = Vec::new();
    for line in text.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 4 || parts[0] != "6R143" || parts[1] != "hap" {
            continue;
        }
        let hash = parts[3].to_string();
        let kv = parse_kv_fields(&parts[4..]);
        out.push(JavaHap {
            hash,
            is_ref: kv.get("isRef").map(|s| s.as_str()) == Some("true"),
            cigar: kv.get("cigar").cloned().unwrap_or_else(|| ".".to_string()),
            loc: kv.get("loc").cloned().unwrap_or_default(),
            align_start: kv.get("alignStart").cloned().unwrap_or_default(),
            bases: kv.get("bases").cloned().unwrap_or_default(),
        });
    }
    out
}

fn load_java_stage_reads(path: &Path, stage: &str) -> BTreeSet<ReadKey> {
    let text = std::fs::read_to_string(path).expect("java dump");
    let mut out = BTreeSet::new();
    for line in text.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 4 || parts[0] != "6R143" || parts[1] != stage {
            continue;
        }
        let qname = parts[2].to_string();
        let kv = parse_kv_fields(&parts[3..]);
        let flags: u16 = kv.get("flags").and_then(|s| s.parse().ok()).unwrap_or(0);
        out.insert((qname, flags));
    }
    out
}

fn load_java_matrix(path: &Path, stage: &str) -> HashMap<CellKey, f64> {
    let text = std::fs::read_to_string(path).expect("java dump");
    let mut out = HashMap::new();
    for line in text.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 4 || parts[0] != "6R143" || parts[1] != stage {
            continue;
        }
        let qname = parts[2].to_string();
        let kv = parse_kv_fields(&parts[3..]);
        let flags: u16 = kv.get("flags").and_then(|s| s.parse().ok()).unwrap_or(0);
        let hap = kv.get("hap").cloned().unwrap_or_default();
        let bits =
            u64::from_str_radix(kv.get("bits").map(|s| s.as_str()).unwrap_or("0"), 16).unwrap_or(0);
        out.insert((qname, flags, hap), f64::from_bits(bits));
    }
    out
}

fn hap_overlap_equal(a: &str, b: &str, read_start: i64, read_end: i64) -> bool {
    if a.is_empty() || b.is_empty() {
        return false;
    }
    let ovl_s = (read_start as u64).max(HAP_START);
    let ovl_e = (read_end as u64).min(HAP_END);
    if ovl_e < ovl_s {
        return true;
    }
    let i0 = (ovl_s - HAP_START) as usize;
    let i1 = (ovl_e - HAP_START) as usize + 1;
    let asl = a.get(i0.min(a.len())..i1.min(a.len())).unwrap_or("");
    let bsl = b.get(i0.min(b.len())..i1.min(b.len())).unwrap_or("");
    asl == bsl
}

#[test]
fn holdout_6r148_best_haplotype_identity() {
    if std::env::var("HOLDOUT_6R148").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R148=1");
        return;
    }
    let _policy = EnvGuard::set("GATK_RS_EXPERIMENTAL_KBEST_POLICY", "unbounded_diagnostic");
    kv(
        "parity_mode",
        "unbounded_diagnostic\tnote=diagnostic_not_production_default",
    );
    kv("resource_policy", "quarantined");
    kv("likelihood_refresh", "NONE");

    let root = repo_root();
    let jf = root.join(JAVA_FLOAT_DUMP_REL);
    let java_haps = load_java_haps(&jf);
    let java_keep = load_java_stage_reads(&jf, "stored");
    let java_stored = load_java_matrix(&jf, "stored");
    assert_eq!(java_haps.len(), 25);
    assert_eq!(java_keep.len(), 201);
    let hap_by_hash: HashMap<String, &JavaHap> =
        java_haps.iter().map(|h| (h.hash.clone(), h)).collect();

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
    let covering = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("ActiveFull")
        .clone();
    let mut java_bounds = covering;
    java_bounds.start = GenomePosition::new_1based(JAVA_ACTIVE_START);
    java_bounds.end = GenomePosition::new_1based(JAVA_ACTIVE_END);
    java_bounds.extended_start = GenomePosition::new_1based(JAVA_PAD_START);
    java_bounds.extended_end = GenomePosition::new_1based(JAVA_PAD_END);

    let args = CallRegionArgs::strict_java();
    begin_likelihood_pipeline_observe();
    begin_poorly_modeled_observe();
    begin_realign_observe();
    let _ = HaplotypeCallerEngine::call_region(&java_bounds, &dict, &ref_fasta, &args)
        .expect("call")
        .expect("outcome");
    let pipe = take_likelihood_pipeline_snaps();
    let _cells = take_likelihood_pipeline_cells();
    let pm = take_poorly_modeled_observe();
    let realign = take_realign_observe();

    let n_pairhmm_snaps = pipe
        .iter()
        .filter(|s| s.stage.contains("pairhmm") || s.stage.contains("refresh"))
        .count();
    kv(
        "post_realign_pairhmm_snaps",
        format!("{n_pairhmm_snaps} (expect 0 named refresh)"),
    );

    let pass0 = pm.iter().map(|r| r.pass).min().unwrap_or(1);
    let rust_keep: BTreeSet<ReadKey> = pm
        .iter()
        .filter(|r| r.pass == pass0 && r.rust_keep)
        .map(|r| (r.qname.clone(), r.flags))
        .collect();
    kv(
        "retained_identity",
        format!(
            "java=201\trust_keep={}\tCOMMON={}\tJAVA_ONLY={}\tRUST_ONLY={}",
            rust_keep.len(),
            java_keep.intersection(&rust_keep).count(),
            java_keep.difference(&rust_keep).count(),
            rust_keep.difference(&java_keep).count()
        ),
    );
    assert_eq!(java_keep.symmetric_difference(&rust_keep).count(), 0);

    let hap_order: Vec<String> = java_haps.iter().map(|h| h.hash.clone()).collect();
    let priorities: Vec<f64> = java_haps
        .iter()
        .map(|h| java_hap_priority(h.is_ref, n_cigar_elements(&h.cigar)))
        .collect();
    let mut java_best: HashMap<ReadKey, (usize, String, f64, Vec<usize>)> = HashMap::new();
    for rk in &java_keep {
        let mut row = vec![f64::NEG_INFINITY; hap_order.len()];
        for (i, h) in hap_order.iter().enumerate() {
            if let Some(&v) = java_stored.get(&(rk.0.clone(), rk.1, h.clone())) {
                row[i] = v;
            }
        }
        let bi = java_search_best_allele(&row, &priorities);
        let best_v = row[bi];
        let ties: Vec<usize> = row
            .iter()
            .enumerate()
            .filter(|(_, &v)| v == best_v)
            .map(|(i, _)| i)
            .collect();
        java_best.insert(rk.clone(), (bi, hap_order[bi].clone(), best_v, ties));
    }

    let rust_best: HashMap<ReadKey, (usize, String)> = realign
        .iter()
        .map(|r| {
            (
                (r.qname.clone(), r.flags),
                (r.best_hap_index, format!("{:016x}", r.best_hap_fnv)),
            )
        })
        .collect();

    let mut mismatches = Vec::new();
    for rk in &java_keep {
        let j = java_best.get(rk).expect("java best");
        let r = rust_best.get(rk).expect("rust best");
        if j.1 != r.1 {
            mismatches.push(rk.clone());
        }
    }
    kv(
        "best_haplotype_selection",
        format!(
            "equal={}/201\tmismatch={}",
            201 - mismatches.len(),
            mismatches.len()
        ),
    );
    assert_eq!(mismatches.len(), 2);

    for rk in &mismatches {
        let (ji, jh, jv, ties) = java_best.get(rk).unwrap();
        let (ri, rh) = rust_best.get(rk).unwrap();
        let jhap = hap_by_hash.get(jh).expect("java hap");
        let rhap = hap_by_hash.get(rh).expect("rust hap");
        let keep_row = realign
            .iter()
            .find(|r| r.qname == rk.0 && r.flags == rk.1)
            .expect("realign row");
        let bases_eq = jhap.bases == rhap.bases;
        let rep_eq = jhap.cigar == rhap.cigar
            && jhap.loc == rhap.loc
            && jhap.align_start == rhap.align_start
            && jhap.is_ref == rhap.is_ref;
        let ovl = hap_overlap_equal(
            &jhap.bases,
            &rhap.bases,
            keep_row.orig_start_1based,
            keep_row.orig_end_1based,
        );
        kv(
            "divergent_read",
            format!(
                "qname={}\tflags={}\tjava_idx={ji}\tjava_hap={jh}\trust_idx={ri}\trust_hap={rh}\tjava_max={jv:.16e}\tjava_tie_n={}\tjava_tie_idx={}\thap_bases_equal={bases_eq}\thap_repr_equal={rep_eq}\toverlap_bases_equal={ovl}\tjava_loc={}\trust_loc={}\tjava_cigar={}\trust_cigar={}",
                rk.0,
                rk.1,
                ties.len(),
                ties.iter().map(|i| i.to_string()).collect::<Vec<_>>().join(","),
                jhap.loc,
                rhap.loc,
                jhap.cigar,
                rhap.cigar
            ),
        );
        assert!(!bases_eq, "6R.147 mismatches are distinct sequences");
        assert!(rep_eq, "same 325M loc/alignStart/isRef");
        assert!(ovl, "overlapping hap span must be identical");
        assert!(!jhap.is_ref && !rhap.is_ref);
    }

    let keep_realign: Vec<_> = realign
        .iter()
        .filter(|r| rust_keep.contains(&(r.qname.clone(), r.flags)))
        .collect();
    assert_eq!(keep_realign.len(), 201);

    let mut bases_eq = 0usize;
    let mut bq_eq = 0usize;
    let mut iq_eq = 0usize;
    let mut dq_eq = 0usize;
    let mut cigar_eq = 0usize;
    let mut start_eq = 0usize;
    let mut end_eq = 0usize;
    let mut flags_eq = 0usize;
    let mut mq_eq = 0usize;
    let mut hc_after = 0usize;
    for r in &keep_realign {
        if r.orig_seq_fnv == r.new_seq_fnv && r.orig_seq_len == r.new_seq_len {
            bases_eq += 1;
        }
        if r.orig_qual_fnv == r.new_qual_fnv {
            bq_eq += 1;
        }
        if r.orig_has_bi == r.new_has_bi && r.orig_bi_fnv == r.new_bi_fnv {
            iq_eq += 1;
        }
        if r.orig_has_bd == r.new_has_bd && r.orig_bd_fnv == r.new_bd_fnv {
            dq_eq += 1;
        }
        if r.orig_cigar == r.new_cigar {
            cigar_eq += 1;
        }
        if r.orig_start_1based == r.new_start_1based {
            start_eq += 1;
        }
        if r.orig_end_1based == r.new_end_1based {
            end_eq += 1;
        }
        if r.flags == r.new_flags {
            flags_eq += 1;
        }
        if r.orig_mq == r.new_mq {
            mq_eq += 1;
        }
        if r.new_has_hc {
            hc_after += 1;
        }
    }
    kv(
        "keep_201_semantic",
        format!(
            "bases={bases_eq}/201\tBQ={bq_eq}/201\tIQ={iq_eq}/201\tDQ={dq_eq}/201\tCIGAR={cigar_eq}/201\tstart={start_eq}/201\tend={end_eq}/201\tflags={flags_eq}/201\tMQ={mq_eq}/201\tHC_after={hc_after}"
        ),
    );
    for rk in &mismatches {
        let r = keep_realign
            .iter()
            .find(|x| x.qname == rk.0 && x.flags == rk.1)
            .expect("mismatch row");
        kv(
            "mismatch_read_state",
            format!(
                "qname={}\tflags={}\tbases={}\tBQ={}\tIQ={}\tDQ={}\tCIGAR={}\tstart={}\tend={}\tflags_eq={}\tMQ={}\tHC_before={}\tHC_after={}",
                rk.0,
                rk.1,
                r.orig_seq_fnv == r.new_seq_fnv,
                r.orig_qual_fnv == r.new_qual_fnv,
                r.orig_bi_fnv == r.new_bi_fnv,
                r.orig_bd_fnv == r.new_bd_fnv,
                r.orig_cigar == r.new_cigar,
                r.orig_start_1based == r.new_start_1based,
                r.orig_end_1based == r.new_end_1based,
                r.flags == r.new_flags,
                r.orig_mq == r.new_mq,
                r.orig_has_hc,
                r.new_has_hc
            ),
        );
        assert_eq!(r.orig_seq_fnv, r.new_seq_fnv);
        assert_eq!(r.orig_qual_fnv, r.new_qual_fnv);
        assert_eq!(r.orig_cigar, r.new_cigar);
        assert_eq!(r.orig_start_1based, r.new_start_1based);
        assert_eq!(r.orig_end_1based, r.new_end_1based);
        assert!(!r.new_has_hc);
    }

    assert_eq!(bases_eq, 201);
    assert_eq!(bq_eq, 201);
    assert_eq!(iq_eq, 201);
    assert_eq!(dq_eq, 201);
    assert_eq!(cigar_eq, 201);
    assert_eq!(start_eq, 201);
    assert_eq!(end_eq, 201);
    assert_eq!(flags_eq, 201);
    assert_eq!(mq_eq, 201);

    kv("classification", "BEST_HAPLOTYPE_IDENTITY_NON_CAUSAL");
    kv(
        "first_divergence",
        "NONE at read/evidence state; winner index leftover is f32/f64 residual",
    );
    kv("production_change", "NONE");
    kv(
        "next_arrow",
        "do not open genotype/AD/PL/QUAL/VCF; winner identity is non-causal at this boundary",
    );
}
