//! 6R.146: `filterPoorlyModeledEvidence` after equivalent 25×249 normalized LLs.
//! Diagnostic-only `unbounded_diagnostic`. Production k-best policy unchanged.
//! Skipped unless `HOLDOUT_6R146=1`.
//!
//! ```text
//! HOLDOUT_6R146=1 GATK_RS_EXPERIMENTAL_KBEST_POLICY=unbounded_diagnostic \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r146_filter_poorly_modeled -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    begin_likelihood_pipeline_observe, begin_poorly_modeled_observe, call_disposition,
    flatten_assembly_regions, take_likelihood_pipeline_cells, take_likelihood_pipeline_snaps,
    take_poorly_modeled_cells, take_poorly_modeled_haplotypes, take_poorly_modeled_observe,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
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
const ESTABLISHED_J_FLOAT_RESIDUAL: f64 = 2.776e-6;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R146\t{key}\t{}", value.as_ref());
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

fn java_log10_min_true_likelihood(qualified_read_len: usize) -> f64 {
    let max_errors = (qualified_read_len as f64 * 0.02).ceil().min(2.0);
    max_errors * -4.0
}

fn java_predicate_keep(max_ll: f64, qlen: usize) -> bool {
    !(max_ll < java_log10_min_true_likelihood(qlen))
}

fn load_java_orig(path: &Path, prefix: &str) -> HashMap<ReadKey, usize> {
    let text = std::fs::read_to_string(path).expect("java dump");
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

fn load_java_stage_reads(path: &Path, prefix: &str, stage: &str) -> BTreeSet<ReadKey> {
    let text = std::fs::read_to_string(path).expect("java dump");
    let mut out = BTreeSet::new();
    for line in text.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 4 || parts[0] != prefix || parts[1] != stage {
            continue;
        }
        let qname = parts[2].to_string();
        let kv = parse_kv_fields(&parts[3..]);
        let flags: u16 = kv.get("flags").and_then(|s| s.parse().ok()).unwrap_or(0);
        out.insert((qname, flags));
    }
    out
}

fn load_java_norm_max(path: &Path, prefix: &str) -> HashMap<ReadKey, (String, f64)> {
    let text = std::fs::read_to_string(path).expect("java dump");
    let mut best: HashMap<ReadKey, (String, f64)> = HashMap::new();
    for line in text.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 4 || parts[0] != prefix || parts[1] != "norm" {
            continue;
        }
        let qname = parts[2].to_string();
        let kv = parse_kv_fields(&parts[3..]);
        let flags: u16 = kv.get("flags").and_then(|s| s.parse().ok()).unwrap_or(0);
        let hap = kv.get("hap").cloned().unwrap_or_default();
        let bits =
            u64::from_str_radix(kv.get("bits").map(|s| s.as_str()).unwrap_or("0"), 16).unwrap_or(0);
        let v = f64::from_bits(bits);
        let rk = (qname, flags);
        match best.get(&rk) {
            Some((_, b)) if v > *b => {
                best.insert(rk, (hap, v));
            }
            None => {
                best.insert(rk, (hap, v));
            }
            _ => {}
        }
    }
    best
}

fn load_java_norm_haps(path: &Path, prefix: &str) -> BTreeSet<String> {
    let text = std::fs::read_to_string(path).expect("java dump");
    let mut out = BTreeSet::new();
    for line in text.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 4 || parts[0] != prefix || parts[1] != "norm" {
            continue;
        }
        let kv = parse_kv_fields(&parts[3..]);
        if let Some(h) = kv.get("hap") {
            out.insert(h.clone());
        }
    }
    out
}

#[test]
fn holdout_6r146_filter_poorly_modeled() {
    if std::env::var("HOLDOUT_6R146").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R146=1");
        return;
    }
    let _policy = EnvGuard::set("GATK_RS_EXPERIMENTAL_KBEST_POLICY", "unbounded_diagnostic");
    kv(
        "parity_mode",
        "unbounded_diagnostic\tnote=diagnostic_not_production_default",
    );
    kv("resource_policy", "quarantined");

    let root = repo_root();
    let jf = root.join(JAVA_FLOAT_DUMP_REL);
    let jd = root.join(JAVA_DOUBLE_DUMP_REL);
    let orig = load_java_orig(&jf, "6R143");
    let java_keep_set = load_java_stage_reads(&jf, "6R143", "stored");
    let java_keep_d = load_java_stage_reads(&jd, "6R144D", "stored");
    let java_norm_reads = load_java_stage_reads(&jf, "6R143", "norm");
    let java_haps = load_java_norm_haps(&jf, "6R143");
    let j_max = load_java_norm_max(&jf, "6R143");
    let d_max = load_java_norm_max(&jd, "6R144D");
    let java_drop: BTreeSet<ReadKey> = orig
        .keys()
        .cloned()
        .filter(|k| !java_keep_set.contains(k))
        .collect();
    kv(
        "java_input",
        format!(
            "orig={}\tnorm_reads={}\tnorm_haps={}\tkeep={}\tdrop={}\tfloat_double_keep_symdiff={}",
            orig.len(),
            java_norm_reads.len(),
            java_haps.len(),
            java_keep_set.len(),
            java_drop.len(),
            java_keep_set.symmetric_difference(&java_keep_d).count()
        ),
    );
    assert_eq!(orig.len(), 249);
    assert_eq!(java_norm_reads.len(), 249);
    assert_eq!(java_haps.len(), 25);
    assert_eq!(java_keep_set.len(), 201);
    assert_eq!(java_drop.len(), 48);
    assert_eq!(java_keep_set, java_keep_d);

    let mut recon_mismatch = 0usize;
    for (rk, &qlen) in &orig {
        let max_ll = j_max.get(rk).map(|(_, v)| *v).unwrap_or(f64::NEG_INFINITY);
        let pred = java_predicate_keep(max_ll, qlen.max(1));
        let live = java_keep_set.contains(rk);
        if pred != live {
            recon_mismatch += 1;
            if recon_mismatch <= 4 {
                kv(
                    "java_recon_mismatch",
                    format!(
                        "qname={}\tflags={}\tqlen={qlen}\tmax_ll={:.16e}\tthr={:.1}\tpred={pred}\tlive={live}",
                        rk.0,
                        rk.1,
                        max_ll,
                        java_log10_min_true_likelihood(qlen.max(1))
                    ),
                );
            }
        }
    }
    kv("java_predicate_vs_stored", recon_mismatch.to_string());
    assert_eq!(
        recon_mismatch, 0,
        "Java stored keep/drop must reconstruct from norm max_ll and orig len"
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
    begin_poorly_modeled_observe();
    let _ = HaplotypeCallerEngine::call_region(&java_bounds, &dict, &ref_fasta, &args)
        .expect("call")
        .expect("outcome");
    let pipe = take_likelihood_pipeline_snaps();
    let cells = take_likelihood_pipeline_cells();
    let pm = take_poorly_modeled_observe();
    let pm_haps = take_poorly_modeled_haplotypes();
    let pm_cells = take_poorly_modeled_cells();
    for s in &pipe {
        kv(
            "pipe_snap",
            format!(
                "seq={}\tstage={}\tn_reads={}\tn_haps={}\tn_ll={}",
                s.seq, s.stage, s.n_reads, s.n_haps, s.n_ll_entries
            ),
        );
    }

    let pass0 = pm.iter().map(|r| r.pass).min().unwrap_or(1);
    let rows: Vec<_> = pm.iter().filter(|r| r.pass == pass0).collect();
    let hap_pass = pm_haps.iter().filter(|h| h.pass == pass0).count();
    let filter_cells = pm_cells.iter().filter(|c| c.pass == pass0).count();
    kv(
        "rust_filter_input",
        format!(
            "rows={}\thap_columns={hap_pass}\tfilter_cells={filter_cells}\tn_hap_cells_min={}\tn_hap_cells_max={}\tn_columns_min={}\textra_retain={}",
            rows.len(),
            rows.iter().map(|r| r.n_hap_cells).min().unwrap_or(0),
            rows.iter().map(|r| r.n_hap_cells).max().unwrap_or(0),
            rows.iter().map(|r| r.n_columns).min().unwrap_or(0),
            rows.iter().filter(|r| r.extra_retain).count()
        ),
    );
    assert_eq!(rows.len(), 249);
    assert_eq!(hap_pass, 25);
    assert!(rows.iter().all(|r| r.n_hap_cells == 25));
    assert!(rows.iter().all(|r| r.n_columns == 25));

    let rust_keys: BTreeSet<ReadKey> = rows.iter().map(|r| (r.qname.clone(), r.flags)).collect();
    assert_eq!(
        rust_keys.len(),
        249,
        "filter rows must be unique QNAME+flags"
    );
    let common = java_norm_reads.intersection(&rust_keys).count();
    let java_only = java_norm_reads.difference(&rust_keys).count();
    let rust_only = rust_keys.difference(&java_norm_reads).count();
    kv(
        "lifecycle_identity",
        format!("COMMON={common}\tJAVA_ONLY={java_only}\tRUST_ONLY={rust_only}"),
    );
    assert_eq!(common, 249);
    assert_eq!(java_only, 0);
    assert_eq!(rust_only, 0);

    let rust_haps: BTreeSet<String> = pm_haps
        .iter()
        .filter(|h| h.pass == pass0)
        .map(|h| format!("{:016x}", h.fnv1a))
        .collect();
    kv(
        "haplotype_population",
        format!(
            "java={}/25\trust={}/25\tCOMMON={}\tJAVA_ONLY={}\tRUST_ONLY={}",
            java_haps.len(),
            rust_haps.len(),
            java_haps.intersection(&rust_haps).count(),
            java_haps.difference(&rust_haps).count(),
            rust_haps.difference(&java_haps).count()
        ),
    );
    assert_eq!(java_haps, rust_haps);

    let seq_norm = cells
        .iter()
        .filter(|c| c.stage == "normalize")
        .map(|c| c.seq)
        .min()
        .expect("normalize");
    let mut norm_by_idx: BTreeMap<usize, ReadKey> = BTreeMap::new();
    for c in cells
        .iter()
        .filter(|c| c.stage == "normalize" && c.seq == seq_norm)
    {
        norm_by_idx
            .entry(c.read_index)
            .or_insert((c.qname.clone(), c.flags));
    }
    let mut lifecycle_idx_mismatch = 0usize;
    for r in &rows {
        let Some(k) = norm_by_idx.get(&r.row_index) else {
            lifecycle_idx_mismatch += 1;
            continue;
        };
        if k.0 != r.qname || k.1 != r.flags {
            lifecycle_idx_mismatch += 1;
        }
    }
    kv(
        "lifecycle_row_i_scored_read_i",
        format!("mismatches={lifecycle_idx_mismatch}/249"),
    );
    assert_eq!(lifecycle_idx_mismatch, 0);

    let mut qlen_eq = 0usize;
    let mut thr_eq = 0usize;
    let mut winner_eq = 0usize;
    let mut winner_in_max_set = 0usize;
    let mut java_pred_eq_rust_java_equiv = 0usize;
    let mut keep_eq = 0usize;
    let mut extra_retain = 0usize;
    let mut could_cross = 0usize;
    let mut max_ll_abs = 0.0f64;
    let mut n_short_thr = 0usize;
    let rust_keep: BTreeSet<ReadKey> = rows
        .iter()
        .filter(|r| r.rust_keep)
        .map(|r| (r.qname.clone(), r.flags))
        .collect();
    for r in &rows {
        let rk = (r.qname.clone(), r.flags);
        let jq = orig.get(&rk).copied().unwrap_or(0).max(1);
        let jthr = java_log10_min_true_likelihood(jq);
        if r.qual_len == jq {
            qlen_eq += 1;
        }
        if (r.threshold - jthr).abs() < 1e-12 {
            thr_eq += 1;
        }
        if jthr > -8.0 {
            n_short_thr += 1;
        }
        let (jwin, jmax) = j_max.get(&rk).cloned().unwrap_or((String::new(), f64::NAN));
        let dmax = d_max.get(&rk).map(|p| p.1).unwrap_or(f64::NAN);
        let rwin = format!("{:016x}", r.argmax_fnv);
        if rwin == jwin {
            winner_eq += 1;
        }
        if (jmax - r.max_ll).abs() <= 1e-12 || (jmax - r.max_ll).abs() < 1.0 {
            // membership in float-max set checked via bits of Java max vs rust hap's Java cell not needed:
            if (jmax - r.max_ll).abs() <= ESTABLISHED_J_FLOAT_RESIDUAL * 2.0 {
                winner_in_max_set += 1;
            }
        }
        let jdec = java_keep_set.contains(&rk);
        if r.java_equiv_keep == jdec {
            java_pred_eq_rust_java_equiv += 1;
        }
        if r.rust_keep == jdec {
            keep_eq += 1;
        } else if keep_eq + extra_retain < 8 {
            kv(
                "keep_mismatch",
                format!(
                    "qname={}\tflags={}\tjqlen={jq}\trqlen={}\tjthr={jthr:.1}\trthr={:.1}\tjmax={:.16e}\trmax={:.16e}\tdmax={:.16e}\tjkeep={jdec}\tr_java_equiv={}\tr_keep={}\textra={}",
                    r.qname,
                    r.flags,
                    r.qual_len,
                    r.threshold,
                    jmax,
                    r.max_ll,
                    dmax,
                    r.java_equiv_keep,
                    r.rust_keep,
                    r.extra_retain
                ),
            );
        }
        if r.extra_retain {
            extra_retain += 1;
        }
        let margin = (jmax - jthr).abs();
        let rmargin = (r.max_ll - r.threshold).abs();
        if margin <= ESTABLISHED_J_FLOAT_RESIDUAL || rmargin <= ESTABLISHED_J_FLOAT_RESIDUAL {
            could_cross += 1;
        }
        let d = (jmax - r.max_ll).abs();
        if d > max_ll_abs {
            max_ll_abs = d;
        }
    }
    kv(
        "qualified_len_threshold",
        format!(
            "qlen_integer_equal={qlen_eq}/249\tthreshold_equal={thr_eq}/249\tjava_threshold_gt_neg8={n_short_thr}"
        ),
    );
    kv(
        "max_ll",
        format!(
            "strict_winner_equal={winner_eq}/249\tmax_ll_within_2x_backend={winner_in_max_set}/249\tmax_abs_delta={max_ll_abs:.16e}\tcould_cross_threshold={could_cross}"
        ),
    );
    kv(
        "keep_drop",
        format!(
            "java_keep={}\trust_keep={}\tequal={keep_eq}/249\textra_retain={extra_retain}\tjava_equiv_pred_equal={java_pred_eq_rust_java_equiv}/249\tsymdiff={}",
            java_keep_set.len(),
            rust_keep.len(),
            java_keep_set.symmetric_difference(&rust_keep).count()
        ),
    );

    let java_keep_order: Vec<ReadKey> = {
        let text = std::fs::read_to_string(&jf).expect("java dump");
        let mut seen = BTreeSet::new();
        let mut order = Vec::new();
        for line in text.lines() {
            let parts: Vec<&str> = line.split('\t').collect();
            if parts.len() < 4 || parts[0] != "6R143" || parts[1] != "stored" {
                continue;
            }
            let qname = parts[2].to_string();
            let kvmap = parse_kv_fields(&parts[3..]);
            let flags: u16 = kvmap.get("flags").and_then(|s| s.parse().ok()).unwrap_or(0);
            let rk = (qname, flags);
            if seen.insert(rk.clone()) {
                order.push(rk);
            }
        }
        order
    };
    let mut rust_keep_order: Vec<(usize, ReadKey)> = rows
        .iter()
        .filter(|r| r.rust_keep)
        .map(|r| (r.row_index, (r.qname.clone(), r.flags)))
        .collect();
    rust_keep_order.sort_by_key(|(i, _)| *i);
    let rust_keep_order: Vec<ReadKey> = rust_keep_order.into_iter().map(|(_, k)| k).collect();
    kv(
        "keep_order",
        format!(
            "java_n={}\trust_n={}\tsequence_equal={}",
            java_keep_order.len(),
            rust_keep_order.len(),
            java_keep_order == rust_keep_order
        ),
    );

    let mut tsv = String::from(
        "QNAME\tflags\tqualifiedLen_Java\tqualifiedLen_Rust\tthreshold_Java\tthreshold_Rust\tmax_ll_J_float\tmax_ll_J_double\tmax_ll_R_f64\twinner_Java\twinner_Rust\tdecision_Java\tdecision_Rust\textra_retain\n",
    );
    let mut sorted_rows = rows.clone();
    sorted_rows.sort_by(|a, b| (&a.qname, a.flags).cmp(&(&b.qname, b.flags)));
    for r in &sorted_rows {
        let rk = (r.qname.clone(), r.flags);
        let jq = orig.get(&rk).copied().unwrap_or(0).max(1);
        let jthr = java_log10_min_true_likelihood(jq);
        let (jwin, jmax) = j_max.get(&rk).cloned().unwrap_or((String::new(), f64::NAN));
        let dmax = d_max.get(&rk).map(|p| p.1).unwrap_or(f64::NAN);
        tsv.push_str(&format!(
            "{}\t{}\t{}\t{}\t{:.1}\t{:.1}\t{:.16e}\t{:.16e}\t{:.16e}\t{}\t{:016x}\t{}\t{}\t{}\n",
            r.qname,
            r.flags,
            jq,
            r.qual_len,
            jthr,
            r.threshold,
            jmax,
            dmax,
            r.max_ll,
            jwin,
            r.argmax_fnv,
            if java_keep_set.contains(&rk) {
                "KEEP"
            } else {
                "DROP"
            },
            if r.rust_keep { "KEEP" } else { "DROP" },
            r.extra_retain
        ));
    }
    let tsv_path = root.join("docs/parity/6R.146_FILTER_ORACLE.tsv");
    std::fs::write(&tsv_path, tsv).expect("write 6R.146 oracle tsv");
    kv("oracle_tsv", tsv_path.display().to_string());
    kv("filter_cells_expected", "6225");
    assert_eq!(filter_cells, 6225);

    assert_eq!(
        thr_eq, 249,
        "threshold must match even if qlen integers differ above 51"
    );
    assert_eq!(java_pred_eq_rust_java_equiv, 249);
    assert_eq!(keep_eq, 249);
    assert_eq!(extra_retain, 0);
    assert_eq!(could_cross, 0);
    assert_eq!(java_keep_set.symmetric_difference(&rust_keep).count(), 0);

    let classification = if java_only != 0 || rust_only != 0 || lifecycle_idx_mismatch != 0 {
        "POORLY_MODELED_EVIDENCE_LIFECYCLE_DIVERGENCE"
    } else if java_haps != rust_haps {
        "POORLY_MODELED_HAPLOTYPE_POPULATION_DIVERGENCE"
    } else if thr_eq != 249 {
        "POORLY_MODELED_THRESHOLD_DIVERGENCE"
    } else if keep_eq != 249 {
        "POORLY_MODELED_THRESHOLD_DIVERGENCE"
    } else {
        "EXPECTED_PRECISION_NON_CAUSAL"
    };
    kv("classification", classification);
    kv(
        "first_divergence",
        if classification == "EXPECTED_PRECISION_NON_CAUSAL" {
            "NONE_at_filterPoorlyModeledEvidence"
        } else {
            "filterPoorlyModeledEvidence"
        },
    );
    kv("production_change", "NONE");
    kv(
        "next_arrow",
        "filterAlleles / realignReadsToTheirBestHaplotype on the retained evidence",
    );
    assert_eq!(classification, "EXPECTED_PRECISION_NON_CAUSAL");
}
