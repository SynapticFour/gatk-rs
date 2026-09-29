//! 6R.294: force only the final trimmed padded span to Java's
//! `20:29455560-29455728` and measure the live genotyping PL.
//! `trim_modern` is not modified. The override is default-off.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r294_trim_span_causal_counterfactual -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    begin_hap_list_observe, begin_likelihood_pipeline_observe, begin_trim_calc_observe,
    call_disposition, flatten_assembly_regions, set_forensic_6r294_padded_span,
    take_colocated_merge_numerics, take_hap_list_snaps, take_hap_list_trim_span,
    take_likelihood_pipeline_cells, take_trim_calc_snap, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, HapListSnap,
    HaplotypeCallerEngine, LikelihoodPipelineCell, ReadFilterParams, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const LOCUS: u64 = 29_455_649;
const ALT: &str = "TGTTTG";
const QNAME: &str = "HISEQ1:11:H8GV6ADXX:1:2116:18670:99941";
const FLAGS: u16 = 99;
const JAVA_SPAN: (u64, u64) = (29_455_560, 29_455_728);
const JAVA_CONTINUOUS: f64 = 3517.46536813194416027;
const CHECK_HAPS: [usize; 5] = [0, 4, 9, 19, 20];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R294\t{key}\t{}", value.as_ref());
}

fn span(bounds: (u64, u64)) -> String {
    format!("{}-{}", bounds.0, bounds.1)
}

fn clear_other_forensics() {
    gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r287_global_hap_floor(None, None);
    gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r288_haplotype_substitute(
        None, None, 0, None,
    );
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_target(None, 0);
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_substitute(None, None);
}

fn continuous_pl(gls: &[f64]) -> Vec<f64> {
    let best = gls.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    gls.iter().map(|g| -10.0 * (g - best)).collect()
}

fn hom_alt_index(allele: usize) -> usize {
    allele * (allele + 1) / 2 + allele
}

fn fmt_f(x: f64) -> String {
    format!("{x:.17}")
}

fn fmt_vec_f(xs: &[f64]) -> String {
    xs.iter().map(|x| fmt_f(*x)).collect::<Vec<_>>().join(",")
}

fn fmt_vec_i(xs: &[i32]) -> String {
    xs.iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn fmt_vec_u(xs: &[u32]) -> String {
    xs.iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

struct Measured {
    production_padded: (u64, u64),
    applied_padded: (u64, u64),
    read_repr: String,
    hap_n: usize,
    hap_fnvs: Vec<u64>,
    hap_geom: String,
    checkpoint: Vec<(usize, f64)>,
    continuous: Vec<f64>,
    hom_alt: f64,
    integer_pl: Vec<i32>,
    integer_hom_alt: i32,
    subset_pl: Vec<i32>,
    emitted_pl: Vec<u32>,
}

fn hap_geom(snaps: &[HapListSnap]) -> (usize, Vec<u64>, String) {
    let snap = snaps
        .iter()
        .find(|s| s.stage == "pairhmm_input")
        .expect("pairhmm_input haplotypes");
    let lens: BTreeSet<_> = snap.columns.iter().map(|c| c.len).collect();
    let spans: BTreeSet<_> = snap
        .columns
        .iter()
        .map(|c| (c.loc_start, c.loc_end))
        .collect();
    let ref_h = snap.columns.iter().find(|c| c.is_reference);
    let geom = match ref_h {
        Some(h) => format!(
            "n={} ref_len={} ref_span={}-{} ref_cigar={} lens={:?} spans={:?}",
            snap.n, h.len, h.loc_start, h.loc_end, h.cigar, lens, spans
        ),
        None => format!("n={} lens={:?} spans={:?}", snap.n, lens, spans),
    };
    let fnvs = snap.columns.iter().map(|c| c.fnv1a).collect();
    (snap.n, fnvs, geom)
}

fn checkpoint(cells: &[LikelihoodPipelineCell]) -> Vec<(usize, f64)> {
    let mut out = Vec::new();
    for hap in CHECK_HAPS {
        let cell = cells.iter().find(|c| {
            c.stage == "post_kernel" && c.qname == QNAME && c.flags == FLAGS && c.hap_index == hap
        });
        if let Some(cell) = cell {
            out.push((hap, cell.log10_likelihood));
        }
    }
    out
}

fn measure(span_override: Option<(u64, u64)>) -> Measured {
    clear_other_forensics();
    set_forensic_6r294_padded_span(span_override.map(|s| s.0), span_override.map(|s| s.1));
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_target(Some(QNAME), FLAGS);
    begin_trim_calc_observe();
    begin_hap_list_observe();
    begin_likelihood_pipeline_observe();

    let root = repo_root();
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
    let region = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= LOCUS
                && r.end.get() >= LOCUS
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

    let trim = take_trim_calc_snap().expect("trim snap");
    let applied = take_hap_list_trim_span().expect("applied span");
    let snaps = take_hap_list_snaps();
    let cells = take_likelihood_pipeline_cells();
    let read_snaps = gatk_haplotypecaller::likelihood_engine::take_forensic_6r289_snaps();
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_target(None, 0);
    let numerics = take_colocated_merge_numerics();
    let site = numerics
        .iter()
        .find(|n| n.loc == LOCUS && n.alts.iter().any(|a| a == ALT))
        .expect("merge numerics");
    let allele = site.alts.iter().position(|a| a == ALT).expect("TGTTTG") + 1;
    let idx = hom_alt_index(allele);
    let continuous = continuous_pl(&site.merged_gls);
    let hom_alt = continuous[idx];
    let integer_pl = site.merged_pl.clone();
    let integer_hom_alt = integer_pl[idx];
    let recs =
        try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let emitted = recs
        .iter()
        .find(|r| r.position == LOCUS && r.alternate.iter().any(|a| a == ALT))
        .expect("emit");
    let read = read_snaps.first().expect("diagnostic read snap");
    let (hap_n, hap_fnvs, hap_geom) = hap_geom(&snaps);
    let measured = Measured {
        production_padded: (trim.padded_start, trim.padded_end),
        applied_padded: (applied.trim_start, applied.trim_end),
        read_repr: format!(
            "{} @ {} len={}",
            read.cigar,
            read.pos_1based,
            read.read_bases.len()
        ),
        hap_n,
        hap_fnvs,
        hap_geom,
        checkpoint: checkpoint(&cells),
        continuous,
        hom_alt,
        integer_pl,
        integer_hom_alt,
        subset_pl: site.subset_pl.clone(),
        emitted_pl: emitted.samples[0].pl.clone().unwrap_or_default(),
    };
    set_forensic_6r294_padded_span(None, None);
    measured
}

#[test]
fn forensic_6r294_trim_span_causal_counterfactual() {
    let root = repo_root();
    if !root.join(REF_REL).is_file() || !root.join(BAM_REL).is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }
    let baseline = measure(None);
    let cf = measure(Some(JAVA_SPAN));

    assert_eq!(baseline.production_padded, (29_455_569, 29_455_724));
    assert_eq!(baseline.applied_padded, baseline.production_padded);
    assert_eq!(cf.production_padded, baseline.production_padded);
    assert_eq!(cf.applied_padded, JAVA_SPAN);

    let hap_same = baseline.hap_fnvs == cf.hap_fnvs && baseline.hap_n == cf.hap_n;
    let gap = (baseline.hom_alt - JAVA_CONTINUOUS).abs();
    let remain = (cf.hom_alt - JAVA_CONTINUOUS).abs();
    let moved = baseline.hom_alt - cf.hom_alt;
    assert!(
        (baseline.hom_alt - 3517.516170057229).abs() < 1e-9,
        "baseline hom-alt {}",
        baseline.hom_alt
    );
    assert!(
        (cf.hom_alt - 3517.46534871093808761).abs() < 1e-9,
        "cf hom-alt {}",
        cf.hom_alt
    );
    assert_eq!(baseline.integer_hom_alt, 3518);
    assert_eq!(cf.integer_hom_alt, 3517);
    assert_eq!(baseline.emitted_pl, vec![570, 0, 3518]);
    assert_eq!(cf.emitted_pl, vec![570, 0, 3517]);
    assert_eq!(baseline.read_repr, "60H88M @ 29455569 len=88");
    assert_eq!(cf.read_repr, "51H97M @ 29455560 len=97");
    let classification = if (baseline.hom_alt - cf.hom_alt).abs() < 1e-3 {
        "TRIM_SPAN_NOT_CAUSAL_FOR_REMAINING_PL"
    } else if remain + 1e-3 < gap {
        "TRIM_SPAN_CAUSAL_FOR_REMAINING_PL"
    } else if remain > gap + 1e-3 {
        "TRIM_SPAN_MOVES_PL_AWAY_FROM_JAVA"
    } else {
        "TRIM_SPAN_NOT_CAUSAL_FOR_REMAINING_PL"
    };

    kv("baseline_trim_span", span(baseline.applied_padded));
    kv("counterfactual_trim_span", span(cf.applied_padded));
    kv("java_trim_span", span(JAVA_SPAN));
    kv("production_trim_unchanged", span(cf.production_padded));
    kv("read_baseline", &baseline.read_repr);
    kv("read_counterfactual", &cf.read_repr);
    kv("read_java", "51H97M @ 29455560 len=97");
    kv("hap_baseline", &baseline.hap_geom);
    kv("hap_counterfactual", &cf.hap_geom);
    kv(
        "hap_java",
        "174bp span 29455560-29455728 (169 reference bases + 5I)",
    );
    kv("hap_set_unchanged", hap_same.to_string());
    kv(
        "checkpoint_baseline",
        baseline
            .checkpoint
            .iter()
            .map(|(h, v)| format!("H{h}={}", fmt_f(*v)))
            .collect::<Vec<_>>()
            .join(" "),
    );
    kv(
        "checkpoint_counterfactual",
        cf.checkpoint
            .iter()
            .map(|(h, v)| format!("H{h}={}", fmt_f(*v)))
            .collect::<Vec<_>>()
            .join(" "),
    );
    kv("baseline_continuous_pl", fmt_vec_f(&baseline.continuous));
    kv("counterfactual_continuous_pl", fmt_vec_f(&cf.continuous));
    kv("baseline_continuous_hom_alt", fmt_f(baseline.hom_alt));
    kv("counterfactual_continuous_hom_alt", fmt_f(cf.hom_alt));
    kv("java_continuous_hom_alt", fmt_f(JAVA_CONTINUOUS));
    kv("baseline_integer_pl", fmt_vec_i(&baseline.integer_pl));
    kv("counterfactual_integer_pl", fmt_vec_i(&cf.integer_pl));
    kv(
        "baseline_integer_hom_alt",
        baseline.integer_hom_alt.to_string(),
    );
    kv(
        "counterfactual_integer_hom_alt",
        cf.integer_hom_alt.to_string(),
    );
    kv("java_integer_hom_alt", "3517");
    kv("baseline_subset_pl", fmt_vec_i(&baseline.subset_pl));
    kv("counterfactual_subset_pl", fmt_vec_i(&cf.subset_pl));
    kv("baseline_emitted_pl", fmt_vec_u(&baseline.emitted_pl));
    kv("counterfactual_emitted_pl", fmt_vec_u(&cf.emitted_pl));
    kv("java_emitted_pl", "570,0,3517");
    kv("continuous_movement_toward_java", fmt_f(gap - remain));
    kv("continuous_delta_baseline_minus_cf", fmt_f(moved));
    assert_eq!(classification, "TRIM_SPAN_CAUSAL_FOR_REMAINING_PL");
    kv("classification", classification);
    kv("production_change", "NONE");
}
