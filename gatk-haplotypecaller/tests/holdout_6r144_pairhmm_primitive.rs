//! 6R.144: primitive PairHMM likelihoods after equivalent input tuples.
//! Diagnostic-only `unbounded_diagnostic`. Production k-best policy unchanged.
//! Skipped unless `HOLDOUT_6R144=1`.
//!
//! Compares identity-aligned (QNAME, flags, haplotype FNV) cells. Does not pin
//! coordinate-specific floating-point bits. Existing 6R.97–6R.99 contracts cover
//! the precision/backend residual.
//!
//! ```text
//! HOLDOUT_6R144=1 GATK_RS_EXPERIMENTAL_KBEST_POLICY=unbounded_diagnostic \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r144_pairhmm_primitive -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    begin_likelihood_pipeline_observe, call_disposition, flatten_assembly_regions,
    resolve_pair_hmm_impl, take_likelihood_pipeline_cells, take_likelihood_pipeline_snaps,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, PairHmmImpl, ReadFilterParams, WalkerTraversalConfig,
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
/// 6R.98/6R.125 established J_double vs Rust scale. Materially larger → DP justified.
const ESTABLISHED_J_DOUBLE_VS_RUST_SCALE: f64 = 1e-8;
const MATERIAL_J_DOUBLE_VS_RUST: f64 = 1e-6;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R144\t{key}\t{}", value.as_ref());
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

fn is_exact_f32_widened(x: f64) -> bool {
    f64::from(x as f32).to_bits() == x.to_bits()
}

fn ordered_f64_bits(x: f64) -> i64 {
    let bits = x.to_bits() as i64;
    if bits < 0 {
        i64::MIN - bits
    } else {
        bits
    }
}

fn ulp_distance(a: f64, b: f64) -> u64 {
    ordered_f64_bits(a).abs_diff(ordered_f64_bits(b))
}

fn java_log10_min_true_likelihood(qualified_read_len: usize) -> f64 {
    let max_errors = (qualified_read_len as f64 * 0.02).ceil().min(2.0);
    max_errors * -4.0
}

fn java_keep(max_ll: f64, qlen: usize) -> bool {
    !(max_ll < java_log10_min_true_likelihood(qlen))
}

type CellKey = (String, u16, String);

#[derive(Clone, Debug)]
struct PrimCell {
    qname: String,
    flags: u16,
    hap: String,
    bits: u64,
    f32wide: bool,
    finite: bool,
}

impl PrimCell {
    fn ll(&self) -> f64 {
        f64::from_bits(self.bits)
    }
}

fn load_java_meta(path: &Path, prefix: &str) -> BTreeMap<String, String> {
    let text = std::fs::read_to_string(path).unwrap_or_default();
    let mut out = BTreeMap::new();
    for line in text.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 || parts[0] != prefix {
            continue;
        }
        if parts[1] == "prim" {
            continue;
        }
        if parts[1].contains('=') {
            continue;
        }
        out.insert(parts[1].to_string(), parts[2..].join("\t"));
    }
    out
}

fn load_java_prim(path: &Path, prefix: &str) -> Vec<PrimCell> {
    let text = std::fs::read_to_string(path).expect("java prim dump");
    let mut out = Vec::new();
    for line in text.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 4 || parts[0] != prefix || parts[1] != "prim" {
            continue;
        }
        let qname = parts[2].to_string();
        let kv = parse_kv_fields(&parts[3..]);
        let flags: u16 = kv.get("flags").and_then(|s| s.parse().ok()).unwrap_or(0);
        let hap = kv.get("hap").cloned().unwrap_or_default();
        let bits =
            u64::from_str_radix(kv.get("bits").map(|s| s.as_str()).unwrap_or("0"), 16).unwrap_or(0);
        let f32wide = kv.get("f32wide").map(|s| s == "true").unwrap_or(false);
        let finite = kv.get("finite").map(|s| s == "true").unwrap_or(false);
        out.push(PrimCell {
            qname,
            flags,
            hap,
            bits,
            f32wide,
            finite,
        });
    }
    out
}

fn load_java_orig_len(path: &Path, prefix: &str) -> HashMap<(String, u16), usize> {
    let text = std::fs::read_to_string(path).expect("java orig dump");
    let mut out = HashMap::new();
    for line in text.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 4 || parts[0] != prefix || parts[1] != "orig" {
            continue;
        }
        let kv = parse_kv_fields(&parts[3..]);
        let qname = kv.get("qname").cloned().unwrap_or_default();
        let flags: u16 = kv.get("flags").and_then(|s| s.parse().ok()).unwrap_or(0);
        let len: usize = kv.get("len").and_then(|s| s.parse().ok()).unwrap_or(0);
        out.insert((qname, flags), len);
    }
    out
}

fn cells_to_map(cells: &[PrimCell]) -> HashMap<CellKey, PrimCell> {
    let mut m = HashMap::new();
    for c in cells {
        m.insert((c.qname.clone(), c.flags, c.hap.clone()), c.clone());
    }
    m
}

#[derive(Debug)]
struct PairStats {
    n: usize,
    exact_equal: usize,
    max_abs: f64,
    mean_abs: f64,
    max_rel: f64,
    max_ulp: u64,
    mean_ulp: f64,
    f32_wide_a: usize,
    f32_wide_b: usize,
}

fn pair_stats(a: &HashMap<CellKey, f64>, b: &HashMap<CellKey, f64>, keys: &[CellKey]) -> PairStats {
    let mut n = 0usize;
    let mut eq = 0usize;
    let mut abs_sum = 0.0;
    let mut max_abs = 0.0;
    let mut max_rel = 0.0;
    let mut ulp_sum = 0u128;
    let mut max_ulp = 0u64;
    let mut f32_a = 0usize;
    let mut f32_b = 0usize;
    for k in keys {
        let Some(&av) = a.get(k) else {
            continue;
        };
        let Some(&bv) = b.get(k) else {
            continue;
        };
        n += 1;
        if is_exact_f32_widened(av) {
            f32_a += 1;
        }
        if is_exact_f32_widened(bv) {
            f32_b += 1;
        }
        if av.to_bits() == bv.to_bits() {
            eq += 1;
        }
        let d = (av - bv).abs();
        abs_sum += d;
        if d > max_abs {
            max_abs = d;
        }
        let denom = av.abs().max(bv.abs()).max(1e-300);
        let rel = d / denom;
        if rel > max_rel {
            max_rel = rel;
        }
        let u = ulp_distance(av, bv);
        ulp_sum += u as u128;
        if u > max_ulp {
            max_ulp = u;
        }
    }
    PairStats {
        n,
        exact_equal: eq,
        max_abs,
        mean_abs: if n == 0 { 0.0 } else { abs_sum / n as f64 },
        max_rel,
        max_ulp,
        mean_ulp: if n == 0 {
            0.0
        } else {
            ulp_sum as f64 / n as f64
        },
        f32_wide_a: f32_a,
        f32_wide_b: f32_b,
    }
}

fn print_stats(label: &str, s: &PairStats) {
    kv(
        "pair",
        format!(
            "label={label}\tn={}\texact_equal={}\tmax_abs={:.16e}\tmean_abs={:.16e}\tmax_rel={:.16e}\tmax_ulp={}\tmean_ulp={:.6}\tf32_wide_a={}\tf32_wide_b={}",
            s.n,
            s.exact_equal,
            s.max_abs,
            s.mean_abs,
            s.max_rel,
            s.max_ulp,
            s.mean_ulp,
            s.f32_wide_a,
            s.f32_wide_b
        ),
    );
}

fn per_read_argmax(
    vals: &HashMap<CellKey, f64>,
    hap_order: &[String],
) -> HashMap<(String, u16), (String, f64, f64)> {
    let mut by_read: BTreeMap<(String, u16), Vec<(String, f64)>> = BTreeMap::new();
    for ((q, f, h), &v) in vals {
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
        let mut second = f64::NEG_INFINITY;
        for (h, v) in cells {
            if v > best_v {
                second = best_v;
                best_v = v;
                best_h = h;
            } else if v > second {
                second = v;
            }
        }
        let margin = if second.is_finite() {
            best_v - second
        } else {
            f64::INFINITY
        };
        out.insert(rk, (best_h, best_v, margin));
    }
    out
}

#[test]
fn holdout_6r144_pairhmm_primitive() {
    if std::env::var("HOLDOUT_6R144").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R144=1");
        return;
    }
    let _policy = EnvGuard::set("GATK_RS_EXPERIMENTAL_KBEST_POLICY", "unbounded_diagnostic");
    kv(
        "parity_mode",
        "unbounded_diagnostic\tnote=diagnostic_not_production_default",
    );
    kv("resource_policy", "quarantined");

    let root = repo_root();
    let java_float_path = root.join(JAVA_FLOAT_DUMP_REL);
    let java_double_path = root.join(JAVA_DOUBLE_DUMP_REL);
    assert!(
        java_float_path.is_file(),
        "missing Java float dump {}",
        java_float_path.display()
    );
    assert!(
        java_double_path.is_file(),
        "missing Java double dump {}",
        java_double_path.display()
    );

    let j_meta = load_java_meta(&java_float_path, "6R143");
    let d_meta = load_java_meta(&java_double_path, "6R144D");
    kv(
        "java_float_path",
        format!(
            "pairhmm_class={}\tkernel_buffer_n={}\tkernel_buffer_f32_wide={}\tmatrix_f32_wide={}\tbuffer_matrix_sorted_bits_equal={}\tprim_n_hap={}\tprim_n_ev={}",
            j_meta.get("pairhmm_class").map(|s| s.as_str()).unwrap_or("."),
            j_meta.get("kernel_buffer_n").map(|s| s.as_str()).unwrap_or("."),
            j_meta.get("kernel_buffer_f32_wide").map(|s| s.as_str()).unwrap_or("."),
            j_meta.get("matrix_f32_wide").map(|s| s.as_str()).unwrap_or("."),
            j_meta.get("buffer_matrix_sorted_bits_equal").map(|s| s.as_str()).unwrap_or("."),
            j_meta.get("prim_n_hap").map(|s| s.as_str()).unwrap_or("."),
            j_meta.get("prim_n_ev").map(|s| s.as_str()).unwrap_or(".")
        ),
    );
    kv(
        "java_double_path",
        format!(
            "pairhmm_class={}\tkernel_buffer_n={}\tkernel_buffer_f32_wide={}\tmatrix_f32_wide={}\tbuffer_matrix_sorted_bits_equal={}\tprim_n_hap={}\tprim_n_ev={}",
            d_meta.get("pairhmm_class").map(|s| s.as_str()).unwrap_or("."),
            d_meta.get("kernel_buffer_n").map(|s| s.as_str()).unwrap_or("."),
            d_meta.get("kernel_buffer_f32_wide").map(|s| s.as_str()).unwrap_or("."),
            d_meta.get("matrix_f32_wide").map(|s| s.as_str()).unwrap_or("."),
            d_meta.get("buffer_matrix_sorted_bits_equal").map(|s| s.as_str()).unwrap_or("."),
            d_meta.get("prim_n_hap").map(|s| s.as_str()).unwrap_or("."),
            d_meta.get("prim_n_ev").map(|s| s.as_str()).unwrap_or(".")
        ),
    );

    let j_cells = load_java_prim(&java_float_path, "6R143");
    let d_cells = load_java_prim(&java_double_path, "6R144D");
    let orig_len = load_java_orig_len(&java_float_path, "6R143");
    kv("java_float_cells", j_cells.len().to_string());
    kv("java_double_cells", d_cells.len().to_string());
    assert_eq!(j_cells.len(), 6225, "Java float prim 25×249");
    assert_eq!(d_cells.len(), 6225, "Java double prim 25×249");

    let rust_backend = resolve_pair_hmm_impl(PairHmmImpl::FastestAvailable);
    kv(
        "rust_backend",
        format!(
            "impl=FastestAvailable\tresolved={}\tpath=score_read_against_haplotypes->score_read_haps_logless->Vec<f64>",
            rust_backend.label()
        ),
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

    let seq0 = cells
        .iter()
        .filter(|c| c.stage == "post_kernel")
        .map(|c| c.seq)
        .min()
        .expect("post_kernel");
    let rust_cells: Vec<_> = cells
        .iter()
        .filter(|c| c.stage == "post_kernel" && c.seq == seq0)
        .collect();
    kv("rust_post_kernel_cells", rust_cells.len().to_string());
    kv(
        "rust_f32_wide",
        rust_cells
            .iter()
            .filter(|c| is_exact_f32_widened(c.log10_likelihood))
            .count()
            .to_string(),
    );

    let j_map = cells_to_map(&j_cells);
    let d_map = cells_to_map(&d_cells);
    let mut r_map: HashMap<CellKey, PrimCell> = HashMap::new();
    for c in &rust_cells {
        r_map.insert(
            (c.qname.clone(), c.flags, format!("{:016x}", c.hap_fnv)),
            PrimCell {
                qname: c.qname.clone(),
                flags: c.flags,
                hap: format!("{:016x}", c.hap_fnv),
                bits: c.log10_likelihood.to_bits(),
                f32wide: is_exact_f32_widened(c.log10_likelihood),
                finite: c.log10_likelihood.is_finite(),
            },
        );
    }

    let j_keys: BTreeSet<CellKey> = j_map.keys().cloned().collect();
    let d_keys: BTreeSet<CellKey> = d_map.keys().cloned().collect();
    let r_keys: BTreeSet<CellKey> = r_map.keys().cloned().collect();
    let matched: Vec<CellKey> = j_keys.intersection(&r_keys).cloned().collect();
    let java_only = j_keys.difference(&r_keys).count();
    let rust_only = r_keys.difference(&j_keys).count();
    let double_only_vs_float = d_keys.symmetric_difference(&j_keys).count();
    kv(
        "primitive_buffer",
        format!(
            "java_float_n={}\tjava_double_n={}\trust_n={}\tmatched={}\tjava_only={}\trust_only={}\tfloat_double_key_symdiff={}",
            j_keys.len(),
            d_keys.len(),
            r_keys.len(),
            matched.len(),
            java_only,
            rust_only,
            double_only_vs_float
        ),
    );
    assert_eq!(j_keys.len(), 6225);
    assert_eq!(d_keys.len(), 6225);
    assert_eq!(matched.len(), 6225);
    assert_eq!(java_only, 0);
    assert_eq!(rust_only, 0);
    assert_eq!(double_only_vs_float, 0);

    let j_ll: HashMap<CellKey, f64> = j_map.iter().map(|(k, c)| (k.clone(), c.ll())).collect();
    let d_ll: HashMap<CellKey, f64> = d_map.iter().map(|(k, c)| (k.clone(), c.ll())).collect();
    let r_ll: HashMap<CellKey, f64> = r_map.iter().map(|(k, c)| (k.clone(), c.ll())).collect();
    let j_wide: HashMap<CellKey, f64> = j_ll
        .iter()
        .map(|(k, v)| (k.clone(), f64::from(*v as f32)))
        .collect();
    let r_wide: HashMap<CellKey, f64> = r_ll
        .iter()
        .map(|(k, v)| (k.clone(), f64::from(*v as f32)))
        .collect();

    let jf_r = pair_stats(&j_ll, &r_ll, &matched);
    let jd_r = pair_stats(&d_ll, &r_ll, &matched);
    let jf_jd = pair_stats(&j_ll, &d_ll, &matched);
    let jf_jwide = pair_stats(&j_ll, &j_wide, &matched);
    let r_rwide = pair_stats(&r_ll, &r_wide, &matched);
    print_stats("J_float_vs_R_f64", &jf_r);
    print_stats("J_double_vs_R_f64", &jd_r);
    print_stats("J_float_vs_J_double", &jf_jd);
    print_stats("J_float_vs_widen_float_J_float", &jf_jwide);
    print_stats("R_f64_vs_widen_float_R_f64", &r_rwide);

    let j_f32 = j_cells.iter().filter(|c| c.f32wide).count();
    let d_f32 = d_cells.iter().filter(|c| c.f32wide).count();
    let j_finite = j_cells.iter().filter(|c| c.finite).count();
    let r_finite = rust_cells
        .iter()
        .filter(|c| c.log10_likelihood.is_finite())
        .count();
    kv(
        "representation",
        format!(
            "J_float_f32wide={j_f32}/6225\tJ_double_f32wide={d_f32}/6225\tJ_float_eq_widen={}/6225\tR_eq_widen={}/6225\tJ_finite={j_finite}\tR_finite={r_finite}",
            jf_jwide.exact_equal, r_rwide.exact_equal
        ),
    );
    assert_eq!(j_f32, 6225, "Java default store is f32-widened");
    assert_eq!(jf_jwide.exact_equal, 6225);
    assert_eq!(
        r_rwide.exact_equal, 0,
        "Rust NeonF64 is not an f32-widened store"
    );
    assert_eq!(d_f32, 0, "Java double store is native f64, not f32-widen");

    let mut hap_order = Vec::new();
    let mut seen = BTreeSet::new();
    for c in &j_cells {
        if seen.insert(c.hap.clone()) {
            hap_order.push(c.hap.clone());
        }
    }
    kv("hap_order_n", hap_order.len().to_string());
    assert_eq!(hap_order.len(), 25);

    let j_win = per_read_argmax(&j_ll, &hap_order);
    let d_win = per_read_argmax(&d_ll, &hap_order);
    let r_win = per_read_argmax(&r_ll, &hap_order);
    assert_eq!(j_win.len(), 249);
    assert_eq!(r_win.len(), 249);

    let mut winner_jf_r = 0usize;
    let mut winner_jd_r = 0usize;
    let mut winner_jf_jd = 0usize;
    let mut min_margin_jf = f64::INFINITY;
    let mut min_margin_jd = f64::INFINITY;
    let mut n_jf_exact_ties = 0usize;
    let mut rust_in_jf_max_set = 0usize;
    let mut rust_in_jd_near_max = 0usize;
    let mut n_keep_eq_jf_r = 0usize;
    let mut n_keep_eq_jd_r = 0usize;
    let mut n_could_cross_jf = 0usize;
    let mut n_could_cross_jd = 0usize;
    let mut missing_len = 0usize;
    let mut n_mismatch_printed = 0usize;
    for (rk, (jh, jmax, jmargin)) in &j_win {
        let rh = r_win.get(rk).expect("rust winner");
        let dh = d_win.get(rk).expect("double winner");
        if jh == &rh.0 {
            winner_jf_r += 1;
        }
        if dh.0 == rh.0 {
            winner_jd_r += 1;
        }
        if jh == &dh.0 {
            winner_jf_jd += 1;
        }
        if *jmargin < min_margin_jf {
            min_margin_jf = *jmargin;
        }
        if dh.2 < min_margin_jd {
            min_margin_jd = dh.2;
        }
        let mut jf_max_set = 0usize;
        let mut rust_jf_eq_max = false;
        let mut rust_jd_near = false;
        for h in &hap_order {
            let k = (rk.0.clone(), rk.1, h.clone());
            if let Some(&jv) = j_ll.get(&k) {
                if jv.to_bits() == jmax.to_bits() {
                    jf_max_set += 1;
                    if h == &rh.0 {
                        rust_jf_eq_max = true;
                    }
                }
            }
            if let Some(&dv) = d_ll.get(&k) {
                if (dh.1 - dv).abs() <= jd_r.max_abs && h == &rh.0 {
                    rust_jd_near = true;
                }
            }
        }
        if jf_max_set > 1 {
            n_jf_exact_ties += 1;
        }
        if rust_jf_eq_max {
            rust_in_jf_max_set += 1;
        }
        if rust_jd_near {
            rust_in_jd_near_max += 1;
        }
        if jh != &rh.0 && n_mismatch_printed < 6 {
            let jk = (rk.0.clone(), rk.1, jh.clone());
            let rk_hap = (rk.0.clone(), rk.1, rh.0.clone());
            kv(
                "winner_mismatch_jf_r",
                format!(
                    "qname={}\tflags={}\tj_hap={}\tr_hap={}\tj_ll_j={:.16e}\tj_ll_r={:.16e}\tr_ll_j={:.16e}\tr_ll_r={:.16e}\tj_margin={:.16e}\tin_jf_max_set={}",
                    rk.0,
                    rk.1,
                    jh,
                    rh.0,
                    j_ll.get(&jk).copied().unwrap_or(f64::NAN),
                    j_ll.get(&rk_hap).copied().unwrap_or(f64::NAN),
                    r_ll.get(&jk).copied().unwrap_or(f64::NAN),
                    rh.1,
                    jmargin,
                    rust_jf_eq_max
                ),
            );
            n_mismatch_printed += 1;
        }
        if dh.0 != rh.0 {
            let dk = (rk.0.clone(), rk.1, dh.0.clone());
            let rk_hap = (rk.0.clone(), rk.1, rh.0.clone());
            kv(
                "winner_mismatch_jd_r",
                format!(
                    "qname={}\tflags={}\td_hap={}\tr_hap={}\td_ll_d={:.16e}\td_ll_r={:.16e}\tr_ll_d={:.16e}\tr_ll_r={:.16e}\td_margin={:.16e}\tin_jd_near_max={}",
                    rk.0,
                    rk.1,
                    dh.0,
                    rh.0,
                    d_ll.get(&dk).copied().unwrap_or(f64::NAN),
                    d_ll.get(&rk_hap).copied().unwrap_or(f64::NAN),
                    r_ll.get(&dk).copied().unwrap_or(f64::NAN),
                    rh.1,
                    dh.2,
                    rust_jd_near
                ),
            );
        }
        let qlen = orig_len.get(rk).copied().unwrap_or(0);
        if qlen == 0 {
            missing_len += 1;
        }
        let jkeep = java_keep(*jmax, qlen.max(1));
        let rkeep = java_keep(rh.1, qlen.max(1));
        let dkeep = java_keep(dh.1, qlen.max(1));
        if jkeep == rkeep {
            n_keep_eq_jf_r += 1;
        }
        if dkeep == rkeep {
            n_keep_eq_jd_r += 1;
        }
        let thr = java_log10_min_true_likelihood(qlen.max(1));
        if (jmax - thr).abs() <= jf_r.max_abs {
            n_could_cross_jf += 1;
        }
        if (dh.1 - thr).abs() <= jd_r.max_abs {
            n_could_cross_jd += 1;
        }
    }
    kv(
        "winner",
        format!(
            "reads=249\tJ_float_vs_R_strict={winner_jf_r}/249\tJ_double_vs_R_strict={winner_jd_r}/249\tJ_float_vs_J_double_strict={winner_jf_jd}/249\trust_in_J_float_max_set={rust_in_jf_max_set}/249\trust_in_J_double_near_max={rust_in_jd_near_max}/249\tJ_float_exact_tie_reads={n_jf_exact_ties}\tmin_J_float_margin={min_margin_jf:.16e}\tmin_J_double_margin={min_margin_jd:.16e}"
        ),
    );
    kv(
        "poorly_modeled_predicate",
        format!(
            "DROP_iff_max_ll_lt_threshold\tJ_float_vs_R_keep={n_keep_eq_jf_r}/249\tJ_double_vs_R_keep={n_keep_eq_jd_r}/249\tcould_cross_J_float_residual={n_could_cross_jf}\tcould_cross_J_double_residual={n_could_cross_jd}\tmissing_orig_len={missing_len}"
        ),
    );
    assert_eq!(n_keep_eq_jf_r, 249);
    assert_eq!(n_keep_eq_jd_r, 249);
    assert_eq!(n_could_cross_jf, 0);
    assert_eq!(n_could_cross_jd, 0);
    assert_eq!(missing_len, 0);
    assert_eq!(
        rust_in_jd_near_max, 249,
        "Rust argmax must lie in the Java-double near-max set (backend residual, not a new ranking)"
    );

    let classification = if jf_r.exact_equal == jf_r.n {
        "UNEXPECTED_EXACT_FLOAT_MATCH"
    } else if jd_r.max_abs > MATERIAL_J_DOUBLE_VS_RUST {
        "UNEXPECTED_PRIMITIVE_DIVERGENCE"
    } else if jd_r.max_abs <= 100.0 * ESTABLISHED_J_DOUBLE_VS_RUST_SCALE
        && n_keep_eq_jf_r == 249
        && rust_in_jd_near_max == 249
    {
        "EXPECTED_PAIRHMM_PRECISION_BACKEND_DELTA"
    } else {
        "UNEXPECTED_PRIMITIVE_DIVERGENCE"
    };
    kv("classification", classification);
    kv(
        "dp_investigation",
        if jd_r.max_abs > MATERIAL_J_DOUBLE_VS_RUST {
            "JUSTIFIED"
        } else {
            "NOT_JUSTIFIED"
        },
    );
    kv("production_change", "NONE");
    kv(
        "first_divergence",
        if classification == "EXPECTED_PAIRHMM_PRECISION_BACKEND_DELTA" {
            "NONE_at_boundary_E_beyond_established_float_vs_f64_backend"
        } else {
            "primitive_PairHMM_result"
        },
    );
    kv(
        "next_arrow",
        "post-kernel normalizeLikelihoods over the equivalent 25×249 matrix",
    );
    kv(
        "established_contract",
        format!(
            "6R.125_this_locus_shared3486_J_float_max_abs~2.465e-6_J_double_max_abs~2.078e-8\t6R.98_seq6_J_double_max_abs=1.7548813957546372e-8\tthis_locus_J_float_max_abs={:.16e}\tthis_locus_J_double_max_abs={:.16e}",
            jf_r.max_abs, jd_r.max_abs
        ),
    );

    assert_eq!(classification, "EXPECTED_PAIRHMM_PRECISION_BACKEND_DELTA");
    assert!(
        jd_r.max_abs <= MATERIAL_J_DOUBLE_VS_RUST,
        "J_double vs R_f64 max_abs {:.3e} exceeds established backend scale",
        jd_r.max_abs
    );
    assert!(jf_r.max_abs > jd_r.max_abs);
    assert!(jf_r.exact_equal < jf_r.n);
}
