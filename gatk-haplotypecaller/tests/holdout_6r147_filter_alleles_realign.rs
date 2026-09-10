//! 6R.147: `filterAlleles` / realign after the equivalent 201×25 poorly-modeled object.
//! Diagnostic-only `unbounded_diagnostic`. Production k-best policy unchanged.
//! Skipped unless `HOLDOUT_6R147=1`.
//!
//! ```text
//! HOLDOUT_6R147=1 GATK_RS_EXPERIMENTAL_KBEST_POLICY=unbounded_diagnostic \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r147_filter_alleles_realign -- --nocapture --test-threads=1
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

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R147\t{key}\t{}", value.as_ref());
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

fn rust_stage_map(
    cells: &[gatk_haplotypecaller::LikelihoodPipelineCell],
    stage: &str,
) -> HashMap<CellKey, f64> {
    let seq = cells
        .iter()
        .filter(|c| c.stage == stage)
        .map(|c| c.seq)
        .min()
        .expect(stage);
    let mut out = HashMap::new();
    for c in cells.iter().filter(|c| c.stage == stage && c.seq == seq) {
        out.insert(
            (c.qname.clone(), c.flags, format!("{:016x}", c.hap_fnv)),
            c.log10_likelihood,
        );
    }
    out
}

#[test]
fn holdout_6r147_filter_alleles_realign() {
    if std::env::var("HOLDOUT_6R147").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R147=1");
        return;
    }
    let _policy = EnvGuard::set("GATK_RS_EXPERIMENTAL_KBEST_POLICY", "unbounded_diagnostic");
    kv(
        "parity_mode",
        "unbounded_diagnostic\tnote=diagnostic_not_production_default",
    );
    kv("resource_policy", "quarantined");
    kv("java_filterAlleles_default", "false");
    kv("java_stepwiseFiltering_default", "false");

    let root = repo_root();
    let jf = root.join(JAVA_FLOAT_DUMP_REL);
    let java_haps = load_java_haps(&jf);
    let java_keep = load_java_stage_reads(&jf, "stored");
    let java_norm = load_java_matrix(&jf, "norm");
    let java_stored = load_java_matrix(&jf, "stored");
    kv(
        "java_input",
        format!(
            "haps={}\tstored_reads={}\tstored_cells={}\tnorm_cells={}",
            java_haps.len(),
            java_keep.len(),
            java_stored.len(),
            java_norm.len()
        ),
    );
    assert_eq!(java_haps.len(), 25);
    assert_eq!(java_keep.len(), 201);
    assert_eq!(java_stored.len(), 201 * 25);
    let java_hap_set: BTreeSet<String> = java_haps.iter().map(|h| h.hash.clone()).collect();
    assert_eq!(java_hap_set.len(), 25);
    assert!(java_haps.iter().all(|h| h.cigar == "325M"));

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
    let cells = take_likelihood_pipeline_cells();
    let pm = take_poorly_modeled_observe();
    let realign = take_realign_observe();
    for s in &pipe {
        kv(
            "pipe_snap",
            format!(
                "seq={}\tstage={}\tn_reads={}\tn_haps={}\tn_ll={}",
                s.seq, s.stage, s.n_reads, s.n_haps, s.n_ll_entries
            ),
        );
    }

    let compaction = pipe.iter().find(|s| s.stage == "compaction");
    let normalize = pipe.iter().find(|s| s.stage == "normalize");
    let rust_haps_at_filter = compaction.map(|s| s.n_haps).unwrap_or(0);
    let rust_reads_at_filter = compaction.map(|s| s.n_reads).unwrap_or(0);
    kv(
        "filter_mask",
        format!(
            "java_skip=identity_25\trust_compaction_haps={rust_haps_at_filter}/25\trust_compaction_reads={rust_reads_at_filter}\tmask_equal={}",
            rust_haps_at_filter == 25
        ),
    );
    assert_eq!(rust_haps_at_filter, 25, "6R.94 68-vs-70 must stay closed");
    assert_eq!(normalize.map(|s| s.n_haps).unwrap_or(0), 25);
    assert_eq!(normalize.map(|s| s.n_reads).unwrap_or(0), 249);

    let pass0 = pm.iter().map(|r| r.pass).min().unwrap_or(1);
    let pm_rows: Vec<_> = pm.iter().filter(|r| r.pass == pass0).collect();
    let rust_keep: BTreeSet<ReadKey> = pm_rows
        .iter()
        .filter(|r| r.rust_keep)
        .map(|r| (r.qname.clone(), r.flags))
        .collect();
    let rust_pm_keys: BTreeSet<ReadKey> =
        pm_rows.iter().map(|r| (r.qname.clone(), r.flags)).collect();
    kv(
        "retained_identity",
        format!(
            "java={}\trust_pm_input={}\trust_keep={}\tCOMMON={}\tJAVA_ONLY={}\tRUST_ONLY={}",
            java_keep.len(),
            rust_pm_keys.len(),
            rust_keep.len(),
            java_keep.intersection(&rust_keep).count(),
            java_keep.difference(&rust_keep).count(),
            rust_keep.difference(&java_keep).count()
        ),
    );
    assert_eq!(java_keep.intersection(&rust_keep).count(), 201);
    assert_eq!(java_keep.difference(&rust_keep).count(), 0);
    assert_eq!(rust_keep.difference(&java_keep).count(), 0);

    let rust_norm = rust_stage_map(&cells, "normalize");
    let rust_hap_set: BTreeSet<String> = rust_norm.keys().map(|k| k.2.clone()).collect();
    kv(
        "haplotype_population",
        format!(
            "java={}/25\trust={}/25\tCOMMON={}\tJAVA_ONLY={}\tRUST_ONLY={}",
            java_hap_set.len(),
            rust_hap_set.len(),
            java_hap_set.intersection(&rust_hap_set).count(),
            java_hap_set.difference(&rust_hap_set).count(),
            rust_hap_set.difference(&java_hap_set).count()
        ),
    );
    assert_eq!(java_hap_set, rust_hap_set);

    let mut n_eq = 0usize;
    let mut n_common_cells = 0usize;
    let mut abs_sum = 0.0f64;
    let mut max_abs = 0.0f64;
    for rk in &java_keep {
        for h in &java_hap_set {
            let jk = (rk.0.clone(), rk.1, h.clone());
            let Some(&jv) = java_stored.get(&jk) else {
                continue;
            };
            let Some(&rv) = rust_norm.get(&jk) else {
                continue;
            };
            n_common_cells += 1;
            if jv.to_bits() == rv.to_bits() {
                n_eq += 1;
            }
            let d = (jv - rv).abs();
            abs_sum += d;
            if d > max_abs {
                max_abs = d;
            }
        }
    }
    kv(
        "likelihood_cells_201x25",
        format!(
            "common={n_common_cells}\texact={n_eq}\tmax_abs={max_abs:.16e}\tmean_abs={:.16e}",
            if n_common_cells == 0 {
                0.0
            } else {
                abs_sum / n_common_cells as f64
            }
        ),
    );
    assert_eq!(n_common_cells, 201 * 25);

    let hap_order: Vec<String> = java_haps.iter().map(|h| h.hash.clone()).collect();
    let priorities: Vec<f64> = java_haps
        .iter()
        .map(|h| java_hap_priority(h.is_ref, n_cigar_elements(&h.cigar)))
        .collect();
    let mut java_best: HashMap<ReadKey, String> = HashMap::new();
    for rk in &java_keep {
        let mut row = vec![f64::NEG_INFINITY; hap_order.len()];
        for (i, h) in hap_order.iter().enumerate() {
            if let Some(&v) = java_stored.get(&(rk.0.clone(), rk.1, h.clone())) {
                row[i] = v;
            }
        }
        let bi = java_search_best_allele(&row, &priorities);
        java_best.insert(rk.clone(), hap_order[bi].clone());
    }

    kv(
        "realign_input",
        format!(
            "java_reads=201\tjava_haps=25\trust_realign_rows={}\trust_realign_n_reads={}\trust_realign_n_haps={}",
            realign.len(),
            realign.first().map(|r| r.n_reads).unwrap_or(0),
            realign.first().map(|r| r.n_haps).unwrap_or(0)
        ),
    );
    assert_eq!(realign.first().map(|r| r.n_haps).unwrap_or(0), 25);
    assert_eq!(realign.first().map(|r| r.n_reads).unwrap_or(0), 249);
    assert_eq!(realign.len(), 249);

    let rust_best: HashMap<ReadKey, String> = realign
        .iter()
        .map(|r| {
            (
                (r.qname.clone(), r.flags),
                format!("{:016x}", r.best_hap_fnv),
            )
        })
        .collect();
    let mut best_eq = 0usize;
    let mut best_ne = 0usize;
    let mut tie_rows = 0usize;
    for rk in &java_keep {
        let jh = java_best.get(rk).cloned().unwrap_or_default();
        let rh = rust_best.get(rk).cloned().unwrap_or_default();
        if jh == rh {
            best_eq += 1;
        } else {
            best_ne += 1;
            if best_ne <= 8 {
                kv(
                    "best_hap_mismatch",
                    format!("qname={}\tflags={}\tjava={jh}\trust={rh}", rk.0, rk.1),
                );
            }
        }
        let mut row = vec![f64::NEG_INFINITY; hap_order.len()];
        for (i, h) in hap_order.iter().enumerate() {
            if let Some(&v) = java_stored.get(&(rk.0.clone(), rk.1, h.clone())) {
                row[i] = v;
            }
        }
        let mut best_v = f64::NEG_INFINITY;
        let mut second_v = f64::NEG_INFINITY;
        for &v in &row {
            if v > best_v {
                second_v = best_v;
                best_v = v;
            } else if v > second_v {
                second_v = v;
            }
        }
        if best_v - second_v < INFORMATIVE {
            tie_rows += 1;
        }
    }
    let rust_norm_best: HashMap<ReadKey, String> = {
        let mut m = HashMap::new();
        for rk in &java_keep {
            let mut row = vec![f64::NEG_INFINITY; hap_order.len()];
            for (i, h) in hap_order.iter().enumerate() {
                if let Some(&v) = rust_norm.get(&(rk.0.clone(), rk.1, h.clone())) {
                    row[i] = v;
                }
            }
            let bi = java_search_best_allele(&row, &priorities);
            m.insert(rk.clone(), hap_order[bi].clone());
        }
        m
    };
    let rust_raw = rust_stage_map(&cells, "compaction");
    let rust_raw_best: HashMap<ReadKey, String> = {
        let mut m = HashMap::new();
        for rk in &java_keep {
            let mut row = vec![f64::NEG_INFINITY; hap_order.len()];
            for (i, h) in hap_order.iter().enumerate() {
                if let Some(&v) = rust_raw.get(&(rk.0.clone(), rk.1, h.clone())) {
                    row[i] = v;
                }
            }
            let bi = java_search_best_allele(&row, &priorities);
            m.insert(rk.clone(), hap_order[bi].clone());
        }
        m
    };
    let mut live_vs_java = 0usize;
    let mut norm_vs_java = 0usize;
    let mut raw_vs_java = 0usize;
    for rk in &java_keep {
        let jh = java_best.get(rk);
        if rust_best.get(rk) == jh {
            live_vs_java += 1;
        }
        if rust_norm_best.get(rk) == jh {
            norm_vs_java += 1;
        }
        if rust_raw_best.get(rk) == jh {
            raw_vs_java += 1;
        }
    }
    kv(
        "best_haplotype_selection",
        format!("live_equal={best_eq}/201\tlive_mismatch={best_ne}\ttie_margin_lt_0.2={tie_rows}"),
    );
    kv(
        "best_hap_causality",
        format!(
            "live_realign_vs_java={live_vs_java}/201\trust_norm_search_vs_java={norm_vs_java}/201\trust_raw_search_vs_java={raw_vs_java}/201"
        ),
    );
    for rk in &java_keep {
        let jh = java_best.get(rk).cloned().unwrap_or_default();
        let rh = rust_best.get(rk).cloned().unwrap_or_default();
        if jh == rh {
            continue;
        }
        let jv = java_stored
            .get(&(rk.0.clone(), rk.1, jh.clone()))
            .copied()
            .unwrap_or(f64::NAN);
        let rv_live_on_java = rust_norm
            .get(&(rk.0.clone(), rk.1, jh.clone()))
            .copied()
            .unwrap_or(f64::NAN);
        let rv_live_on_rust = rust_norm
            .get(&(rk.0.clone(), rk.1, rh.clone()))
            .copied()
            .unwrap_or(f64::NAN);
        kv(
            "best_hap_mismatch_ll",
            format!(
                "qname={}\tflags={}\tjava_hap={jh}\trust_hap={rh}\tjava_ll_on_java_hap={jv:.16e}\trust_norm_on_java_hap={rv_live_on_java:.16e}\trust_norm_on_rust_hap={rv_live_on_rust:.16e}\tnorm_search={}\traw_search={}",
                rk.0,
                rk.1,
                rust_norm_best.get(rk).cloned().unwrap_or_default(),
                rust_raw_best.get(rk).cloned().unwrap_or_default()
            ),
        );
    }

    let keep_realign: Vec<_> = realign
        .iter()
        .filter(|r| rust_keep.contains(&(r.qname.clone(), r.flags)))
        .collect();
    let cigar_changed = keep_realign
        .iter()
        .filter(|r| r.orig_cigar != r.new_cigar || r.orig_start_1based != r.new_start_1based)
        .count();
    kv(
        "realignment_keep_reads",
        format!(
            "n={}\tcigar_or_start_changed={cigar_changed}\tunchanged={}",
            keep_realign.len(),
            keep_realign.len() - cigar_changed
        ),
    );

    let classification = if java_keep.symmetric_difference(&rust_keep).count() != 0 {
        "EVIDENCE_COMPACTION_DIVERGENCE"
    } else if rust_haps_at_filter != 25 || java_hap_set != rust_hap_set {
        "FILTER_ALLELES_DIVERGENCE"
    } else if best_ne != 0 {
        "BEST_HAPLOTYPE_SELECTION_DIVERGENCE"
    } else {
        "EXPECTED_PRECISION_NON_CAUSAL"
    };
    kv("classification", classification);
    kv(
        "first_divergence",
        if classification == "EXPECTED_PRECISION_NON_CAUSAL" {
            "NONE_at_filterAlleles; realign best-hap of KEEP 201 matches Java searchBestAllele"
        } else {
            classification
        },
    );
    kv(
        "likelihood_refresh",
        "NONE (Java stepwiseFiltering=false; Rust change_evidence identity)",
    );
    kv("production_change", "NONE");
    kv(
        "next_arrow",
        "do not open marginalize; 2 KEEP reads disagree on best haplotype identity (CIGAR/start of KEEP 201 unchanged)",
    );
    assert_eq!(rust_haps_at_filter, 25);
    assert_eq!(java_keep.symmetric_difference(&rust_keep).count(), 0);
    assert_eq!(classification, "BEST_HAPLOTYPE_SELECTION_DIVERGENCE");
    assert_eq!(best_ne, 2);
    assert_eq!(cigar_changed, 0);
}
