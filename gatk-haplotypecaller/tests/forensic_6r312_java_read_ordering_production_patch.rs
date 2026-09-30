//! 6R.312: production clipped-read order is Java `ReadCoordinateComparator`.
//! The pre-6R.312 `(tid, pos, qname)` sort remains available only when that
//! diagnostic switch is off, for the causality check below.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r312_java_read_ordering_production_patch -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::assembly_region_finalize::{
    clip_finalized_reads_in_place, java_read_coordinate_compare,
    set_forensic_6r308_java_read_coordinate_order,
};
use gatk_haplotypecaller::read_assembly_filter::{passes_assembly_read, AssemblyReadFilterConfig};
use gatk_haplotypecaller::{
    begin_hap_list_observe, begin_trim_calc_observe, call_disposition, failed_mate_key_of,
    flatten_assembly_regions, passes_mate_on_same_contig_or_no_mapped_mate,
    set_forensic_6r294_padded_span, set_forensic_6r306_exclude_failed_mate,
    take_colocated_merge_numerics, take_hap_list_snaps, take_hap_list_trim_span,
    take_likelihood_pipeline_cells, take_trim_calc_snap, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegion, AssemblyRegionCallDisposition, CallRegionArgs,
    GenomePosition, HapListSnap, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use rust_htslib::bam::record::{Cigar, CigarString};
use rust_htslib::bam::{HeaderView, Record};
use std::collections::{BTreeSet, HashMap};
use std::path::{Path, PathBuf};
use std::sync::Arc;

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_TSV: &str = "gatk-haplotypecaller/tests/6r297_java_hap_ll.tsv";
const LOCUS: u64 = 29_455_649;
const ALT: &str = "TGTTTG";
const JAVA_SPAN: (u64, u64) = (29_455_560, 29_455_728);
const DIAG_QNAME: &str = "HISEQ1:11:H8GV6ADXX:1:2116:18670:99941";
const DIAG_FLAGS: u16 = 99;
const MATE_QNAME: &str = "HISEQ1:11:H8GV6ADXX:2:1103:14252:55237";
const MATE_FLAGS: u16 = 97;
const FORWARD_QNAME: &str = "HISEQ1:11:H8GV6ADXX:1:1116:4033:65919";
const FORWARD_FLAGS: u16 = 99;
const REVERSE_QNAME: &str = "HISEQ1:11:H8GV6ADXX:1:1109:9994:7054";
const REVERSE_FLAGS: u16 = 83;
const JAVA_GL: [&str; 6] = [
    "-517.27823987412966744",
    "-475.12025841147004712",
    "-808.78825973442883424",
    "-460.29478556820799895",
    "-736.51068053854339723",
    "-812.04132043930155760",
];
const JAVA_PL: [&str; 6] = [
    "569.83454305921668492",
    "148.25472843262048173",
    "3484.93474166220858024",
    "-0.00000000000000000",
    "2762.15894970335375547",
    "3517.46534871093535912",
];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R312\t{key}\t{}", value.as_ref());
}

fn fmt_f(x: f64) -> String {
    format!("{x:.17}")
}

fn fmt_vec(xs: &[f64]) -> String {
    xs.iter().map(|x| fmt_f(*x)).collect::<Vec<_>>().join(",")
}

fn cell_bits(cell: &str) -> f64 {
    let hex = cell.split("bits=").nth(1).unwrap().trim();
    f64::from_bits(u64::from_str_radix(hex.trim_start_matches("0x"), 16).unwrap())
}

fn field_cells(text: &str, key: &str) -> Vec<f64> {
    for line in text.lines() {
        let mut p = line.split('\t');
        if p.next() == Some("6R297") && p.next() == Some(key) {
            return p.map(cell_bits).collect();
        }
    }
    panic!("missing {key}");
}

fn ulp_distance(a: f64, b: f64) -> u64 {
    if a.to_bits() == b.to_bits() {
        return 0;
    }
    if !a.is_finite() || !b.is_finite() {
        return u64::MAX;
    }
    fn ordered(x: f64) -> u64 {
        let bits = x.to_bits();
        if bits & (1 << 63) != 0 {
            !bits
        } else {
            bits | (1 << 63)
        }
    }
    ordered(a).abs_diff(ordered(b))
}

fn continuous_pl(gls: &[f64]) -> Vec<f64> {
    let best = gls.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    gls.iter().map(|g| -10.0 * (g - best)).collect()
}

fn header() -> Arc<HeaderView> {
    Arc::new(HeaderView::from_bytes(
        b"@HD\tVN:1.0\n@SQ\tSN:chrHold\tLN:2000\n",
    ))
}

fn placed(name: &str, flags: u16, pos_1based: i64) -> Record {
    let mut r = Record::new();
    r.set_header(header());
    r.set(
        name.as_bytes(),
        Some(&CigarString::from(vec![Cigar::Match(30)])),
        b"ACGTACGTACGTACGTACGTACGTACGTAC",
        &vec![30u8; 30],
    );
    r.set_tid(0);
    r.set_pos(pos_1based - 1);
    r.set_flags(flags);
    r.set_mapq(40);
    r
}

fn tiny_region() -> AssemblyRegion {
    AssemblyRegion {
        contig: "chrHold".into(),
        start: GenomePosition::new_1based(400),
        end: GenomePosition::new_1based(600),
        is_active: true,
        extended_start: GenomePosition::new_1based(1),
        extended_end: GenomePosition::new_1based(2000),
        extension: 100,
        reads: Vec::new(),
        read_qnames: Vec::new(),
        reference: gatk_haplotypecaller::ReferenceContext::empty(),
        features: gatk_haplotypecaller::FeatureContext::empty(),
        pileup_loci: Vec::new(),
    }
}

#[test]
fn production_clip_sort_puts_forward_before_reverse() {
    set_forensic_6r308_java_read_coordinate_order(true);
    let forward = placed("hold-fwd", 99, 500);
    let reverse = placed("a-rev", 83, 500);
    assert!(!forward.is_reverse());
    assert!(reverse.is_reverse());
    assert!(java_read_coordinate_compare(&forward, &reverse) < 0);
    assert_eq!(
        java_read_coordinate_compare(&forward, &reverse),
        -java_read_coordinate_compare(&reverse, &forward)
    );
    let mut reads = vec![reverse, forward];
    clip_finalized_reads_in_place(&mut reads, &tiny_region());
    assert_eq!(reads[0].qname(), b"hold-fwd");
    assert_eq!(reads[1].qname(), b"a-rev");
    let again = reads.clone();
    let mut second = reads.clone();
    clip_finalized_reads_in_place(&mut second, &tiny_region());
    assert_eq!(
        second
            .iter()
            .map(|r| r.qname().to_vec())
            .collect::<Vec<_>>(),
        again.iter().map(|r| r.qname().to_vec()).collect::<Vec<_>>()
    );
}

#[test]
fn comparator_tie_breakers_stay_in_java_order() {
    let mut low = placed("SAME", 99, 500);
    let mut high = placed("SAME", 163, 500);
    assert_eq!(java_read_coordinate_compare(&low, &high), -1);
    low.set_flags(99);
    high.set_flags(99);
    low.set_mapq(10);
    high.set_mapq(50);
    assert_eq!(java_read_coordinate_compare(&low, &high), -1);
    low.set_mapq(50);
    high.set_mapq(50);
    low.set_mtid(0);
    high.set_mtid(0);
    low.set_mpos(4);
    high.set_mpos(20);
    assert_eq!(java_read_coordinate_compare(&low, &high), -1);
    high.set_mpos(4);
    low.set_insert_size(8);
    high.set_insert_size(40);
    assert_eq!(java_read_coordinate_compare(&low, &high), -1);
    assert_eq!(java_read_coordinate_compare(&low, &low), 0);
}

#[test]
fn pre_6r312_sort_keeps_name_order_at_the_same_start() {
    set_forensic_6r308_java_read_coordinate_order(false);
    let mut reads = vec![placed("a-rev", 83, 500), placed("hold-fwd", 99, 500)];
    clip_finalized_reads_in_place(&mut reads, &tiny_region());
    assert_eq!(reads[0].qname(), b"a-rev");
    assert_eq!(reads[1].qname(), b"hold-fwd");
    set_forensic_6r308_java_read_coordinate_order(true);
}

struct Pass {
    order: Vec<(String, u16)>,
    alleles: HashMap<(String, u16), [f64; 3]>,
    gl: Vec<f64>,
    pl: Vec<f64>,
    integer: String,
    subset: String,
    mate_cells: usize,
    mate_in_rows: bool,
    span: (u64, u64),
    geometry_ok: bool,
    active: usize,
    passing: usize,
    mate_fail: usize,
}

fn one_pass(root: &Path, java_order: bool) -> Pass {
    gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r287_global_hap_floor(None, None);
    gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r288_haplotype_substitute(
        None, None, 0, None,
    );
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_substitute(None, None);
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_target(
        Some(DIAG_QNAME),
        DIAG_FLAGS,
    );
    set_forensic_6r294_padded_span(None, None);
    set_forensic_6r306_exclude_failed_mate(false);
    set_forensic_6r308_java_read_coordinate_order(java_order);
    begin_trim_calc_observe();
    begin_hap_list_observe();
    gatk_haplotypecaller::begin_likelihood_pipeline_observe();

    let ref_fasta = root.join(REF_REL);
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, INTERVAL).expect("interval");
    let walk = traverse_assembly_region_walker(
        &dict,
        &specs,
        &ref_fasta,
        &root.join(BAM_REL),
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
    let filter_cfg = AssemblyReadFilterConfig::gatk_defaults();
    let active = region.reads.len();
    let passing = region
        .reads
        .iter()
        .filter(|r| passes_assembly_read(r, &filter_cfg))
        .count();
    let mate_fail = region
        .reads
        .iter()
        .filter(|r| failed_mate_key_of(r).is_some())
        .count();
    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let trim = take_trim_calc_snap().expect("trim");
    let applied = take_hap_list_trim_span().expect("applied");
    let snaps = take_hap_list_snaps();
    let cells = take_likelihood_pipeline_cells();
    let read_snaps = gatk_haplotypecaller::likelihood_engine::take_forensic_6r289_snaps();
    let numerics = take_colocated_merge_numerics();
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_target(None, 0);

    let site = numerics
        .iter()
        .find(|n| n.loc == LOCUS && n.alts.iter().any(|a| a == ALT))
        .expect("site");
    let col = |name: &str| -> usize {
        if site.long_ref == name {
            0
        } else {
            site.alts.iter().position(|a| a == name).unwrap() + 1
        }
    };
    let cols = [col("T"), col("TTTG"), col(ALT)];
    let mut order = Vec::new();
    let mut alleles = HashMap::new();
    for i in 0..site.ad_row_qname.len() {
        let key = (site.ad_row_qname[i].clone(), site.ad_row_flags[i]);
        let row = &site.ad_row_lls[i];
        alleles.insert(key.clone(), [row[cols[0]], row[cols[1]], row[cols[2]]]);
        order.push(key);
    }
    let hap = snaps
        .iter()
        .find(|s: &&HapListSnap| s.stage == "pairhmm_input")
        .expect("haps");
    let ref_h = hap.columns.iter().find(|c| c.is_reference).expect("ref");
    let lens: BTreeSet<_> = hap.columns.iter().map(|c| c.len).collect();
    let read = read_snaps.first().expect("diag");
    let read_repr = format!(
        "{} @ {} len={}",
        read.cigar,
        read.pos_1based,
        read.read_bases.len()
    );
    let geometry_ok = (trim.padded_start, trim.padded_end) == JAVA_SPAN
        && (applied.trim_start, applied.trim_end) == JAVA_SPAN
        && read_repr == "51H97M @ 29455560 len=97"
        && ref_h.len == 169
        && ref_h.cigar == "169M"
        && lens == BTreeSet::from([169, 170, 172, 174])
        && hap.n == 24;
    let integer = site
        .merged_pl
        .iter()
        .map(|p| p.to_string())
        .collect::<Vec<_>>()
        .join(",");
    let recs =
        try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let emitted = recs
        .iter()
        .find(|r| r.position == LOCUS && r.alternate.iter().any(|a| a == ALT))
        .expect("emit");
    let subset = emitted.samples[0]
        .pl
        .clone()
        .unwrap_or_default()
        .iter()
        .map(|p| p.to_string())
        .collect::<Vec<_>>()
        .join(",");
    let gl = site.merged_gls.clone();
    Pass {
        mate_cells: cells
            .iter()
            .filter(|c| c.qname == MATE_QNAME && c.flags == MATE_FLAGS)
            .count(),
        mate_in_rows: order.iter().any(|k| k.0 == MATE_QNAME && k.1 == MATE_FLAGS),
        span: (applied.trim_start, applied.trim_end),
        geometry_ok,
        active,
        passing,
        mate_fail,
        pl: continuous_pl(&gl),
        gl,
        integer,
        subset,
        order,
        alleles,
    }
}

fn load_java(text: &str) -> (Vec<(String, u16)>, HashMap<(String, u16), [f64; 3]>) {
    let mut order = Vec::new();
    let mut alleles = HashMap::new();
    for line in text.lines() {
        let p: Vec<_> = line.split('\t').collect();
        if p.len() >= 7 && p[0] == "6R297" && p[1] == "allele_ll" {
            let key = (p[2].to_string(), p[3].parse().unwrap());
            let vals: Vec<f64> = p[4..].iter().copied().map(cell_bits).collect();
            order.push(key.clone());
            alleles.insert(key, [vals[0], vals[1], vals[2]]);
        }
    }
    (order, alleles)
}

#[test]
fn forensic_6r312_java_read_ordering_production_patch() {
    let root = repo_root();
    let text = std::fs::read_to_string(root.join(JAVA_TSV)).unwrap();
    let (java, java_alleles) = load_java(&text);
    let java_gl = field_cells(&text, "java_gl_vector");
    let java_pl = field_cells(&text, "java_continuous_pl_vector");
    assert_eq!(java.len(), 122);
    assert_eq!(java[0].0, FORWARD_QNAME);
    assert_eq!(java[0].1, FORWARD_FLAGS);
    assert_eq!(java[41].0, REVERSE_QNAME);
    assert_eq!(java[41].1, REVERSE_FLAGS);
    for i in 0..6 {
        assert_eq!(fmt_f(java_gl[i]), JAVA_GL[i]);
        assert_eq!(fmt_f(java_pl[i]), JAVA_PL[i]);
    }
    if !root.join(REF_REL).is_file() || !root.join(BAM_REL).is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }

    let old = one_pass(&root, false);
    let production = one_pass(&root, true);
    set_forensic_6r308_java_read_coordinate_order(true);

    let old_ulp: Vec<u64> = old
        .gl
        .iter()
        .zip(java_gl.iter())
        .map(|(r, j)| ulp_distance(*r, *j))
        .collect();
    let prod_ulp: Vec<u64> = production
        .gl
        .iter()
        .zip(java_gl.iter())
        .map(|(r, j)| ulp_distance(*r, *j))
        .collect();
    let positional = production
        .order
        .iter()
        .zip(java.iter())
        .filter(|(r, j)| r != j)
        .count();
    let mut compared = 0usize;
    let mut differing = 0usize;
    let mut max_abs = 0.0_f64;
    let mut max_ulp = 0u64;
    for key in &java {
        let rust = production.alleles[key];
        let java_row = java_alleles[key];
        let old_row = old.alleles[key];
        for a in 0..3 {
            compared += 1;
            let abs = (rust[a] - java_row[a]).abs();
            let ulp = ulp_distance(rust[a], java_row[a]);
            max_abs = max_abs.max(abs);
            if ulp != u64::MAX {
                max_ulp = max_ulp.max(ulp);
            }
            if rust[a].to_bits() != java_row[a].to_bits()
                || rust[a].to_bits() != old_row[a].to_bits()
            {
                differing += 1;
            }
        }
    }
    let gl_match = production
        .gl
        .iter()
        .zip(java_gl.iter())
        .all(|(r, j)| r.to_bits() == j.to_bits());
    let pl_match = production
        .pl
        .iter()
        .zip(java_pl.iter())
        .all(|(r, j)| r.to_bits() == j.to_bits());
    let order_ok = production.order == java && positional == 0 && production.order.len() == 122;
    let causality_ok = old.order[0].0 == REVERSE_QNAME
        && old.order[0].1 == REVERSE_FLAGS
        && old.order != java
        && old_ulp == [0, 0, 0, 2, 0, 1]
        && compared == 366
        && differing == 0
        && max_ulp == 0;
    let membership_ok = production.active == 450
        && production.passing == 449
        && production.mate_fail == 1
        && !production.mate_in_rows
        && production.mate_cells == 0
        && !old.mate_in_rows
        && old.mate_cells == 0;
    let canonical_still_fails = {
        let mut rec = placed(MATE_QNAME, MATE_FLAGS, 29_455_546);
        rec.set_mtid(1);
        rec.set_paired();
        rec.unset_unmapped();
        rec.unset_mate_unmapped();
        !passes_mate_on_same_contig_or_no_mapped_mate(&rec)
    };
    let classification = if order_ok
        && causality_ok
        && membership_ok
        && production.geometry_ok
        && old.geometry_ok
        && gl_match
        && pl_match
        && production.integer == "570,148,3485,0,2762,3517"
        && production.subset == "570,0,3517"
        && canonical_still_fails
    {
        "JAVA_READ_ORDERING_PRODUCTION_PATCH_REPRODUCES_JAVA"
    } else if order_ok && !gl_match {
        "JAVA_READ_ORDERING_PRODUCTION_PATCH_PARTIAL"
    } else if production.order.len() != 122 || production.mate_in_rows {
        "JAVA_READ_ORDERING_PRODUCTION_PATCH_REGRESSION"
    } else {
        "JAVA_READ_ORDERING_PRODUCTION_PATCH_BLOCKED"
    };

    kv("java_reads", java.len().to_string());
    kv("rust_reads", production.order.len().to_string());
    kv(
        "identity_permutation",
        (production.order == java).to_string(),
    );
    kv("positional_differences", positional.to_string());
    kv(
        "index0",
        format!("{} {}", production.order[0].0, production.order[0].1),
    );
    kv(
        "index41",
        format!("{} {}", production.order[41].0, production.order[41].1),
    );
    kv("compared_cells", compared.to_string());
    kv("differing_cells", differing.to_string());
    kv("max_abs", fmt_f(max_abs));
    kv("max_ulp", max_ulp.to_string());
    kv("old_gl_ulp", format!("{old_ulp:?}"));
    kv("production_gl_ulp", format!("{prod_ulp:?}"));
    kv("final_gl", fmt_vec(&production.gl));
    kv("final_pl", fmt_vec(&production.pl));
    kv("integer_pl", &production.integer);
    kv("emitted_subset_pl", &production.subset);
    kv(
        "span",
        format!("{}-{}", production.span.0, production.span.1),
    );
    kv(
        "failed_mate_absent",
        (!production.mate_in_rows && production.mate_cells == 0).to_string(),
    );
    kv("classification", classification);

    assert!(causality_ok, "old-order ULP {old_ulp:?}");
    assert!(membership_ok, "membership");
    assert!(production.geometry_ok && old.geometry_ok, "STR geometry");
    assert!(order_ok, "positional differences {positional}");
    assert!(gl_match, "{}", fmt_vec(&production.gl));
    assert!(pl_match, "{}", fmt_vec(&production.pl));
    assert_eq!(
        classification,
        "JAVA_READ_ORDERING_PRODUCTION_PATCH_REPRODUCES_JAVA"
    );
}
