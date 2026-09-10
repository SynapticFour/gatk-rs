//! 6R.145: post-kernel `normalizeLikelihoods` after equivalent 25×249 primitives.
//! Diagnostic-only `unbounded_diagnostic`. Production k-best policy unchanged.
//! Skipped unless `HOLDOUT_6R145=1`.
//!
//! ```text
//! HOLDOUT_6R145=1 GATK_RS_EXPERIMENTAL_KBEST_POLICY=unbounded_diagnostic \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r145_normalize -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    begin_likelihood_pipeline_observe, call_disposition, flatten_assembly_regions,
    take_likelihood_pipeline_cells, take_likelihood_pipeline_snaps,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, LikelihoodPipelineCell, ReadFilterParams, WalkerTraversalConfig,
};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_FLOAT_DUMP_REL: &str = "parity/giab/runs/local-pairhmm-diff/6r143_java.txt";
const JAVA_DOUBLE_DUMP_REL: &str = "parity/giab/runs/local-pairhmm-diff/6r144_java_double.txt";
const TARGET: u64 = 29_456_196;
const JAVA_ACTIVE_START: u64 = 29_455_995;
const JAVA_ACTIVE_END: u64 = 29_456_293;
const JAVA_PAD_START: u64 = 29_455_895;
const JAVA_PAD_END: u64 = 29_456_393;
const CAP: f64 = -4.5;
const ESTABLISHED_J_DOUBLE_VS_RUST: f64 = 2.078_263_605_653_774_0e-8;
const MATERIAL_NEW_DIVERGENCE: f64 = 1e-3;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R145\t{key}\t{}", value.as_ref());
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

type CellKey = (String, u16, String);

fn load_java_matrix(path: &Path, prefix: &str, stage: &str) -> HashMap<CellKey, f64> {
    let text = std::fs::read_to_string(path).expect("java dump");
    let mut out = HashMap::new();
    for line in text.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 4 || parts[0] != prefix || parts[1] != stage {
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

fn load_java_meta_count(path: &Path, prefix: &str, key: &str) -> usize {
    let text = std::fs::read_to_string(path).unwrap_or_default();
    for line in text.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 3 && parts[0] == prefix && parts[1] == key {
            return parts[2].parse().unwrap_or(0);
        }
    }
    0
}

fn rust_stage_map(cells: &[LikelihoodPipelineCell], stage: &str) -> HashMap<CellKey, f64> {
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

fn per_read_max(m: &HashMap<CellKey, f64>) -> HashMap<(String, u16), (String, f64)> {
    let mut best: HashMap<(String, u16), (String, f64)> = HashMap::new();
    for ((q, f, h), &v) in m {
        let rk = (q.clone(), *f);
        match best.get(&rk) {
            Some((_, bv)) if v > *bv => {
                best.insert(rk, (h.clone(), v));
            }
            None => {
                best.insert(rk, (h.clone(), v));
            }
            Some((bh, bv)) if v == *bv && h.as_str() < bh.as_str() => {
                // keep first by list-order later; hash-order is only a placeholder
                let _ = (bh, bv);
            }
            _ => {}
        }
    }
    best
}

fn per_read_max_list_order(
    m: &HashMap<CellKey, f64>,
    hap_order: &[String],
) -> HashMap<(String, u16), (String, f64)> {
    let mut by_read: BTreeMap<(String, u16), Vec<(String, f64)>> = BTreeMap::new();
    for ((q, f, h), &v) in m {
        by_read
            .entry((q.clone(), *f))
            .or_default()
            .push((h.clone(), v));
    }
    let mut out = HashMap::new();
    for (rk, mut cells) in by_read {
        cells.sort_by(|a, b| {
            let ia = hap_order
                .iter()
                .position(|h| h == &a.0)
                .unwrap_or(usize::MAX);
            let ib = hap_order
                .iter()
                .position(|h| h == &b.0)
                .unwrap_or(usize::MAX);
            ia.cmp(&ib)
        });
        let mut best_h = String::new();
        let mut best_v = f64::NEG_INFINITY;
        for (h, v) in cells {
            if v > best_v {
                best_v = v;
                best_h = h;
            }
        }
        out.insert(rk, (best_h, best_v));
    }
    out
}

/// Java `normalizeLikelihoodsPerEvidence`: `max(raw, max_P + cap)` with cap = −4.5.
fn predict_normalize(prim: &HashMap<CellKey, f64>) -> HashMap<CellKey, f64> {
    let maxes = per_read_max(prim);
    let mut out = HashMap::new();
    for (k, &raw) in prim {
        let rk = (k.0.clone(), k.1);
        let max_p = maxes.get(&rk).map(|(_, v)| *v).unwrap_or(raw);
        let floor = max_p + CAP;
        let v = if raw < floor { floor } else { raw };
        out.insert(k.clone(), v);
    }
    out
}

fn floored(prim: &HashMap<CellKey, f64>) -> HashMap<CellKey, bool> {
    let maxes = per_read_max(prim);
    let mut out = HashMap::new();
    for (k, &raw) in prim {
        let rk = (k.0.clone(), k.1);
        let max_p = maxes.get(&rk).map(|(_, v)| *v).unwrap_or(raw);
        let floor = max_p + CAP;
        out.insert(k.clone(), raw < floor);
    }
    out
}

fn pair_stats(
    a: &HashMap<CellKey, f64>,
    b: &HashMap<CellKey, f64>,
    keys: &[CellKey],
) -> (usize, usize, f64, f64) {
    let mut n = 0usize;
    let mut eq = 0usize;
    let mut abs_sum = 0.0;
    let mut max_abs = 0.0;
    for k in keys {
        let Some(&av) = a.get(k) else {
            continue;
        };
        let Some(&bv) = b.get(k) else {
            continue;
        };
        n += 1;
        if av.to_bits() == bv.to_bits() {
            eq += 1;
        }
        let d = (av - bv).abs();
        abs_sum += d;
        if d > max_abs {
            max_abs = d;
        }
    }
    (
        n,
        eq,
        max_abs,
        if n == 0 { 0.0 } else { abs_sum / n as f64 },
    )
}

fn hap_order_from(m: &HashMap<CellKey, f64>) -> Vec<String> {
    let mut seen = BTreeSet::new();
    let mut order = Vec::new();
    let mut keys: Vec<_> = m.keys().cloned().collect();
    keys.sort();
    for (_, _, h) in keys {
        if seen.insert(h.clone()) {
            order.push(h);
        }
    }
    order
}

#[test]
fn holdout_6r145_normalize() {
    if std::env::var("HOLDOUT_6R145").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R145=1");
        return;
    }
    let _policy = EnvGuard::set("GATK_RS_EXPERIMENTAL_KBEST_POLICY", "unbounded_diagnostic");
    kv(
        "parity_mode",
        "unbounded_diagnostic\tnote=diagnostic_not_production_default",
    );
    kv("resource_policy", "quarantined");
    kv("cap", format!("{CAP}"));

    let root = repo_root();
    let java_float_path = root.join(JAVA_FLOAT_DUMP_REL);
    let java_double_path = root.join(JAVA_DOUBLE_DUMP_REL);

    let j_prim = load_java_matrix(&java_float_path, "6R143", "prim");
    let j_norm = load_java_matrix(&java_float_path, "6R143", "norm");
    let d_prim = load_java_matrix(&java_double_path, "6R144D", "prim");
    let d_norm = load_java_matrix(&java_double_path, "6R144D", "norm");
    let j_norm_n_hap = load_java_meta_count(&java_float_path, "6R143", "norm_n_hap");
    let j_norm_n_ev = load_java_meta_count(&java_float_path, "6R143", "norm_n_ev");
    let d_norm_n_hap = load_java_meta_count(&java_double_path, "6R144D", "norm_n_hap");
    kv(
        "java_dump",
        format!(
            "float_prim={}\tfloat_norm={}\tfloat_norm_n_hap={j_norm_n_hap}\tfloat_norm_n_ev={j_norm_n_ev}\tdouble_prim={}\tdouble_norm={}\tdouble_norm_n_hap={d_norm_n_hap}",
            j_prim.len(),
            j_norm.len(),
            d_prim.len(),
            d_norm.len()
        ),
    );
    assert_eq!(j_prim.len(), 6225);
    assert_eq!(j_norm.len(), 6225);
    assert_eq!(d_prim.len(), 6225);
    assert_eq!(d_norm.len(), 6225);
    assert_eq!(j_norm_n_hap, 25);
    assert_eq!(d_norm_n_hap, 25);
    assert_eq!(j_norm_n_ev, 249);

    let j_pred = predict_normalize(&j_prim);
    let d_pred = predict_normalize(&d_prim);
    let j_keys: Vec<CellKey> = j_prim.keys().cloned().collect();
    let (n_jp, eq_jp, max_jp, _) = pair_stats(&j_pred, &j_norm, &j_keys);
    let d_keys: Vec<CellKey> = d_prim.keys().cloned().collect();
    let (n_dp, eq_dp, max_dp, _) = pair_stats(&d_pred, &d_norm, &d_keys);
    kv(
        "java_formula_reconstruction",
        format!(
            "J_float_pred_vs_dumped_norm n={n_jp} exact={eq_jp} max_abs={max_jp:.16e}\tJ_double_pred_vs_dumped_norm n={n_dp} exact={eq_dp} max_abs={max_dp:.16e}\tformula=max(raw,max_P-4.5)"
        ),
    );
    assert_eq!(
        eq_jp, 6225,
        "Java float dumped norm must be max(raw, max_P-4.5)"
    );
    assert_eq!(
        eq_dp, 6225,
        "Java double dumped norm must be max(raw, max_P-4.5)"
    );

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
    let _ = HaplotypeCallerEngine::call_region(&java_bounds, &dict, &ref_fasta, &args)
        .expect("call")
        .expect("outcome");
    let pipe = take_likelihood_pipeline_snaps();
    let cells = take_likelihood_pipeline_cells();
    for s in &pipe {
        kv(
            "pipe_snap",
            format!(
                "seq={}\tstage={}\tn_reads={}\tn_haps={}\tn_ll={}",
                s.seq, s.stage, s.n_reads, s.n_haps, s.n_ll_entries
            ),
        );
    }
    let rust_prim = rust_stage_map(&cells, "post_kernel");
    let rust_norm = rust_stage_map(&cells, "normalize");
    kv(
        "rust_matrix",
        format!(
            "post_kernel={}\tnormalize={}",
            rust_prim.len(),
            rust_norm.len()
        ),
    );

    let j_haps: BTreeSet<String> = j_norm.keys().map(|k| k.2.clone()).collect();
    let r_haps: BTreeSet<String> = rust_norm.keys().map(|k| k.2.clone()).collect();
    let common_haps = j_haps.intersection(&r_haps).count();
    let java_only_haps = j_haps.difference(&r_haps).count();
    let rust_only_haps = r_haps.difference(&j_haps).count();
    kv(
        "normalization_population",
        format!(
            "java={}/25\trust={}/25\tCOMMON={common_haps}\tJAVA_ONLY={java_only_haps}\tRUST_ONLY={rust_only_haps}",
            j_haps.len(),
            r_haps.len()
        ),
    );
    assert_eq!(j_haps.len(), 25);
    assert_eq!(r_haps.len(), 25);
    assert_eq!(common_haps, 25);
    assert_eq!(java_only_haps, 0);
    assert_eq!(rust_only_haps, 0);

    let rust_pred = predict_normalize(&rust_prim);
    let r_keys: Vec<CellKey> = rust_prim.keys().cloned().collect();
    let (n_rp, eq_rp, max_rp, _) = pair_stats(&rust_pred, &rust_norm, &r_keys);
    kv(
        "rust_formula_reconstruction",
        format!(
            "post_kernel_all25_pred_vs_captured_normalize n={n_rp} exact={eq_rp} max_abs={max_rp:.16e}\tnote=mismatch_would_mean_restricted_P_mask"
        ),
    );
    assert_eq!(
        eq_rp, 6225,
        "Rust captured normalize must equal all-25-column max(raw, max_P-4.5); EventMap mask would diverge"
    );

    let hap_order = hap_order_from(&j_prim);
    let j_max = per_read_max_list_order(&j_prim, &hap_order);
    let d_max = per_read_max_list_order(&d_prim, &hap_order);
    let r_max = per_read_max_list_order(&rust_prim, &hap_order);
    assert_eq!(j_max.len(), 249);
    assert_eq!(r_max.len(), 249);

    let mut winner_jf_r = 0usize;
    let mut winner_jd_r = 0usize;
    let mut rust_in_jf_max_set = 0usize;
    let mut max_p_abs = 0.0f64;
    let mut max_p_jd_abs = 0.0f64;
    for (rk, (jh, jv)) in &j_max {
        let rh = r_max.get(rk).expect("rust max");
        let dh = d_max.get(rk).expect("double max");
        if jh == &rh.0 {
            winner_jf_r += 1;
        }
        if dh.0 == rh.0 {
            winner_jd_r += 1;
        }
        let mut in_set = false;
        for h in &hap_order {
            let k = (rk.0.clone(), rk.1, h.clone());
            if let Some(&jv_h) = j_prim.get(&k) {
                if jv_h.to_bits() == jv.to_bits() && h == &rh.0 {
                    in_set = true;
                }
            }
        }
        if in_set {
            rust_in_jf_max_set += 1;
        }
        let d1 = (jv - rh.1).abs();
        if d1 > max_p_abs {
            max_p_abs = d1;
        }
        let d2 = (dh.1 - rh.1).abs();
        if d2 > max_p_jd_abs {
            max_p_jd_abs = d2;
        }
    }
    kv(
        "max_p",
        format!(
            "reads=249\tJ_float_vs_R_strict_winner={winner_jf_r}/249\tJ_double_vs_R_strict_winner={winner_jd_r}/249\trust_in_J_float_max_set={rust_in_jf_max_set}/249\tmax_abs_delta_J_float_vs_R={max_p_abs:.16e}\tmax_abs_delta_J_double_vs_R={max_p_jd_abs:.16e}"
        ),
    );
    assert_eq!(rust_in_jf_max_set, 249);

    let j_floor = floored(&j_prim);
    let r_floor = floored(&rust_prim);
    let d_floor = floored(&d_prim);
    let matched: Vec<CellKey> = j_prim
        .keys()
        .filter(|k| rust_prim.contains_key(*k))
        .cloned()
        .collect();
    assert_eq!(matched.len(), 6225);
    let mut floor_eq_jf_r = 0usize;
    let mut floor_eq_jd_r = 0usize;
    let mut n_floored_j = 0usize;
    let mut n_floored_r = 0usize;
    for k in &matched {
        let jf = *j_floor.get(k).unwrap();
        let rf = *r_floor.get(k).unwrap();
        let df = *d_floor.get(k).unwrap();
        if jf {
            n_floored_j += 1;
        }
        if rf {
            n_floored_r += 1;
        }
        if jf == rf {
            floor_eq_jf_r += 1;
        }
        if df == rf {
            floor_eq_jd_r += 1;
        }
    }
    kv(
        "floor_decisions",
        format!(
            "predicate=raw_lt_max_P_minus_4.5\tJ_float_vs_R={floor_eq_jf_r}/6225\tJ_double_vs_R={floor_eq_jd_r}/6225\tjava_floored={n_floored_j}\trust_floored={n_floored_r}"
        ),
    );
    assert_eq!(floor_eq_jf_r, 6225);
    assert_eq!(floor_eq_jd_r, 6225);

    let (n_nr, eq_nr, max_nr, mean_nr) = pair_stats(&j_norm, &rust_norm, &matched);
    let (n_dr, eq_dr, max_dr, mean_dr) = pair_stats(&d_norm, &rust_norm, &matched);
    let (n_pr, eq_pr, max_pr, mean_pr) = pair_stats(&j_prim, &rust_prim, &matched);
    kv(
        "normalized_matrix",
        format!(
            "J_float_vs_R n={n_nr} exact={eq_nr} max_abs={max_nr:.16e} mean_abs={mean_nr:.16e}\tJ_double_vs_R n={n_dr} exact={eq_dr} max_abs={max_dr:.16e} mean_abs={mean_dr:.16e}"
        ),
    );
    kv(
        "primitive_residual_unchanged",
        format!(
            "J_float_vs_R n={n_pr} exact={eq_pr} max_abs={max_pr:.16e} mean_abs={mean_pr:.16e}"
        ),
    );
    assert_eq!(n_nr, 6225);
    assert_eq!(
        eq_nr, 0,
        "J_float vs R must not be bit-identical (backend contract)"
    );
    assert!(
        max_nr < MATERIAL_NEW_DIVERGENCE,
        "normalized max_abs {max_nr} looks like a population/formula split, not backend residual"
    );
    assert!(
        max_dr <= 100.0 * ESTABLISHED_J_DOUBLE_VS_RUST,
        "J_double vs R normalized residual {max_dr} exceeds established backend scale"
    );

    let classification = if java_only_haps != 0 || rust_only_haps != 0 || r_haps.len() != 25 {
        "NORMALIZATION_POPULATION_DIVERGENCE"
    } else if eq_rp != 6225 || eq_jp != 6225 {
        "NORMALIZATION_FORMULA_DIVERGENCE"
    } else if floor_eq_jf_r == 6225 && max_nr < MATERIAL_NEW_DIVERGENCE {
        "EXPECTED_NORMALIZATION_PRECISION_DELTA"
    } else {
        "NORMALIZATION_FORMULA_DIVERGENCE"
    };
    kv("classification", classification);
    kv(
        "first_divergence",
        if classification == "EXPECTED_NORMALIZATION_PRECISION_DELTA" {
            "NONE_at_normalize_beyond_established_pairhmm_backend"
        } else {
            "normalizeLikelihoods"
        },
    );
    kv("production_change", "NONE");
    kv(
        "next_arrow",
        "filterPoorlyModeledEvidence on the normalized 25×249 object",
    );
    kv(
        "regression_6r126_6r129",
        "java_P=25 rust_P=25 JAVA_ONLY=0 RUST_ONLY=0 captured_normalize=all25_formula EventMap_mask_not_used",
    );
    assert_eq!(classification, "EXPECTED_NORMALIZATION_PRECISION_DELTA");
}
