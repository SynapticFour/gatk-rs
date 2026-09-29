//! 6R.311: production exclusion of a paired read whose mapped mate is on
//! another contig. The existing mate predicate is unchanged.
//! This test forces the pre-6R.312 clipped-read sort so its order
//! assertion remains that historical measurement.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r311_failed_mate_production_patch -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::assembly_region_finalize::set_forensic_6r308_java_read_coordinate_order;
use gatk_haplotypecaller::hc_genotyping_engine::diploid_genotype_log10_likelihoods_from_allele_rows;
use gatk_haplotypecaller::read_assembly_filter::{passes_assembly_read, AssemblyReadFilterConfig};
use gatk_haplotypecaller::ReadLikelihoodRow;
use gatk_haplotypecaller::{
    begin_hap_list_observe, begin_likelihood_pipeline_observe, begin_trim_calc_observe,
    call_disposition, failed_mate_evidence_keys, failed_mate_key_of, flatten_assembly_regions,
    forensic_6r306_mate_fail_noted, passes_mate_on_same_contig_or_no_mapped_mate,
    set_forensic_6r294_padded_span, set_forensic_6r306_exclude_failed_mate,
    take_colocated_merge_numerics, take_hap_list_snaps, take_hap_list_trim_span,
    take_likelihood_pipeline_cells, take_trim_calc_snap, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, HapListSnap,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
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
const RUST_FIRST: &str = "HISEQ1:11:H8GV6ADXX:1:1109:9994:7054";
const RUST_FIRST_FLAGS: u16 = 83;
const JAVA_FIRST: &str = "HISEQ1:11:H8GV6ADXX:1:1116:4033:65919";
const JAVA_FIRST_FLAGS: u16 = 99;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R311\t{key}\t{}", value.as_ref());
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
        b"@HD\tVN:1.0\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:1000\n",
    ))
}

fn named(qname: &[u8], paired: bool, unmapped: bool, mate_unmapped: bool, mtid: i32) -> Record {
    let mut r = Record::new();
    r.set_header(header());
    r.set(
        qname,
        Some(&CigarString::from(vec![Cigar::Match(10)])),
        b"AAAAAAAAAA",
        &vec![30u8; 10],
    );
    r.set_tid(0);
    r.set_mtid(mtid);
    if paired {
        r.set_paired();
    }
    if unmapped {
        r.set_unmapped();
    } else {
        r.unset_unmapped();
    }
    if mate_unmapped {
        r.set_mate_unmapped();
    } else {
        r.unset_mate_unmapped();
    }
    r
}

#[test]
fn predicate_matches_java_membership_cases() {
    let unpaired = named(b"unpaired", false, false, false, 1);
    let mate_unmapped = named(b"mate-unmapped", true, false, true, 1);
    let read_unmapped = named(b"read-unmapped", true, true, false, 1);
    let same = named(b"same-contig", true, false, false, 0);
    let other = named(b"other-contig", true, false, false, 1);
    assert!(passes_mate_on_same_contig_or_no_mapped_mate(&unpaired));
    assert!(passes_mate_on_same_contig_or_no_mapped_mate(&mate_unmapped));
    assert!(passes_mate_on_same_contig_or_no_mapped_mate(&read_unmapped));
    assert!(passes_mate_on_same_contig_or_no_mapped_mate(&same));
    assert!(!passes_mate_on_same_contig_or_no_mapped_mate(&other));

    let originals = [unpaired, mate_unmapped, read_unmapped, same, other];
    let keys = failed_mate_evidence_keys(&originals);
    assert_eq!(keys.len(), 1);
    assert!(keys.contains(&(b"other-contig".to_vec(), originals[4].flags())));
    let mut evidence = originals.to_vec();
    evidence.retain(|r| !keys.contains(&(r.qname().to_vec(), r.flags())));
    assert_eq!(evidence.len(), 4);
    assert!(evidence.iter().all(|r| r.qname() != b"other-contig"));
    assert!(evidence.iter().any(|r| r.qname() == b"same-contig"));
    assert!(evidence.iter().any(|r| r.qname() == b"mate-unmapped"));
    assert!(evidence.iter().any(|r| r.qname() == b"unpaired"));
    assert!(evidence.iter().any(|r| r.qname() == b"read-unmapped"));
}

struct JavaAllele {
    order: Vec<(String, u16)>,
    alleles: HashMap<(String, u16), [f64; 3]>,
}

fn load_java(text: &str) -> JavaAllele {
    let mut order = Vec::new();
    let mut alleles = HashMap::new();
    for line in text.lines() {
        let p: Vec<_> = line.split('\t').collect();
        if p.len() < 5 || p[0] != "6R297" || p[1] != "allele_ll" {
            continue;
        }
        let qname = p[2].to_string();
        let flags: u16 = p[3].parse().unwrap();
        let vals: Vec<f64> = p[4..].iter().copied().map(cell_bits).collect();
        order.push((qname.clone(), flags));
        alleles.insert((qname, flags), [vals[0], vals[1], vals[2]]);
    }
    JavaAllele { order, alleles }
}

fn pairhmm_snap(snaps: &[HapListSnap]) -> &HapListSnap {
    snaps
        .iter()
        .find(|s| s.stage == "pairhmm_input")
        .expect("pairhmm_input haplotypes")
}

#[test]
fn forensic_6r311_failed_mate_production_patch() {
    let text = std::fs::read_to_string(repo_root().join(JAVA_TSV)).unwrap();
    let java = load_java(&text);
    let java_gl = field_cells(&text, "java_gl_vector");
    assert_eq!(java.order.len(), 122);
    assert_eq!(java.order[0].0, JAVA_FIRST);
    assert_eq!(java.order[0].1, JAVA_FIRST_FLAGS);

    for rel in [
        "src/engine.rs",
        "src/engine_likelihoods.rs",
        "src/engine_observe.rs",
        "src/read_pre_mate.rs",
    ] {
        let src = std::fs::read_to_string(Path::new(env!("CARGO_MANIFEST_DIR")).join(rel)).unwrap();
        assert!(!src.contains(MATE_QNAME), "{rel}");
        assert!(!src.contains("mate_contig == 7"), "{rel}");
    }

    gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r287_global_hap_floor(None, None);
    gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r288_haplotype_substitute(
        None, None, 0, None,
    );
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_target(
        Some(DIAG_QNAME),
        DIAG_FLAGS,
    );
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_substitute(None, None);
    set_forensic_6r294_padded_span(None, None);
    set_forensic_6r306_exclude_failed_mate(false);
    set_forensic_6r308_java_read_coordinate_order(false);
    begin_trim_calc_observe();
    begin_hap_list_observe();
    begin_likelihood_pipeline_observe();

    let root = repo_root();
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
    let active_before = region.reads.len();
    let passing_filter = region
        .reads
        .iter()
        .filter(|r| passes_assembly_read(r, &filter_cfg))
        .count();
    let mate_fail_before = region
        .reads
        .iter()
        .filter(|r| failed_mate_key_of(r).is_some())
        .count();
    let canonical_fails = region.reads.iter().any(|r| {
        std::str::from_utf8(r.qname()).ok() == Some(MATE_QNAME)
            && r.flags() == MATE_FLAGS
            && !passes_mate_on_same_contig_or_no_mapped_mate(r)
    });
    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let mate_fail_noted = forensic_6r306_mate_fail_noted();
    let trim = take_trim_calc_snap().expect("trim");
    let applied = take_hap_list_trim_span().expect("applied");
    let snaps = take_hap_list_snaps();
    let cells = take_likelihood_pipeline_cells();
    let read_snaps = gatk_haplotypecaller::likelihood_engine::take_forensic_6r289_snaps();
    let numerics = take_colocated_merge_numerics();
    set_forensic_6r306_exclude_failed_mate(false);
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_target(None, 0);

    let site = numerics
        .iter()
        .find(|n| n.loc == LOCUS && n.alts.iter().any(|a| a == ALT))
        .expect("site");
    let mate_in_genotyping = outcome.genotyping_reads.iter().any(|r| {
        std::str::from_utf8(r.qname()).ok() == Some(MATE_QNAME) && r.flags() == MATE_FLAGS
    });
    let mate_row = site
        .ad_row_qname
        .iter()
        .zip(site.ad_row_flags.iter())
        .any(|(q, f)| q == MATE_QNAME && *f == MATE_FLAGS);
    let mate_cells = cells
        .iter()
        .filter(|c| c.qname == MATE_QNAME && c.flags == MATE_FLAGS)
        .count();
    let pairhmm_keys: BTreeSet<(String, u16)> = cells
        .iter()
        .filter(|c| c.stage == "post_kernel")
        .map(|c| (c.qname.clone(), c.flags))
        .collect();

    let column_of = |name: &str| -> Option<usize> {
        if site.long_ref == name {
            return Some(0);
        }
        site.alts.iter().position(|a| a == name).map(|i| i + 1)
    };
    let cols = [
        column_of("T").unwrap(),
        column_of("TTTG").unwrap(),
        column_of(ALT).unwrap(),
    ];
    let mut rust_by_key: HashMap<(String, u16), [f64; 3]> = HashMap::new();
    for i in 0..site.ad_row_lls.len() {
        let row = &site.ad_row_lls[i];
        rust_by_key.insert(
            (site.ad_row_qname[i].clone(), site.ad_row_flags[i]),
            [row[cols[0]], row[cols[1]], row[cols[2]]],
        );
    }
    let java_only = java
        .order
        .iter()
        .filter(|k| !rust_by_key.contains_key(*k))
        .count();
    let rust_only = rust_by_key
        .keys()
        .filter(|k| !java.alleles.contains_key(*k))
        .count();
    let mut compared = 0usize;
    let mut differing = 0usize;
    let mut max_abs = 0.0_f64;
    let mut max_ulp = 0u64;
    for key in &java.order {
        let rust = rust_by_key[key];
        let java_row = java.alleles[key];
        for a in 0..3 {
            compared += 1;
            let abs = (rust[a] - java_row[a]).abs();
            let ulp = ulp_distance(rust[a], java_row[a]);
            if abs > max_abs {
                max_abs = abs;
            }
            if ulp != u64::MAX && ulp > max_ulp {
                max_ulp = ulp;
            }
            if rust[a].to_bits() != java_row[a].to_bits() {
                differing += 1;
            }
        }
    }

    let mut java_order_rows = Vec::new();
    for key in &java.order {
        java_order_rows.push(ReadLikelihoodRow {
            read_index: java_order_rows.len(),
            read_id: String::new(),
            haplotype_log10_likelihoods: rust_by_key[key].to_vec(),
        });
    }
    let replay = diploid_genotype_log10_likelihoods_from_allele_rows(&java_order_rows, 3);
    let rust_gl = site.merged_gls.clone();
    let rust_pl = continuous_pl(&rust_gl);
    let gl_ulp: Vec<u64> = rust_gl
        .iter()
        .zip(java_gl.iter())
        .map(|(r, j)| ulp_distance(*r, *j))
        .collect();
    let order_is_old = site.ad_row_qname.first().map(String::as_str) == Some(RUST_FIRST)
        && site.ad_row_flags.first().copied() == Some(RUST_FIRST_FLAGS)
        && site
            .ad_row_qname
            .iter()
            .position(|q| q == JAVA_FIRST)
            .unwrap()
            > site
                .ad_row_qname
                .iter()
                .position(|q| q == RUST_FIRST)
                .unwrap();

    let hap = pairhmm_snap(&snaps);
    let ref_h = hap.columns.iter().find(|c| c.is_reference).unwrap();
    let lens: BTreeSet<_> = hap.columns.iter().map(|c| c.len).collect();
    let read = read_snaps.first().expect("diagnostic read");
    let read_repr = format!(
        "{} @ {} len={}",
        read.cigar,
        read.pos_1based,
        read.read_bases.len()
    );
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
    let emitted_pl = emitted.samples[0]
        .pl
        .clone()
        .unwrap_or_default()
        .iter()
        .map(|p| p.to_string())
        .collect::<Vec<_>>()
        .join(",");

    let membership_ok = active_before == 450
        && passing_filter == 449
        && mate_fail_before == 1
        && mate_fail_noted == 1
        && canonical_fails
        && !mate_in_genotyping
        && !mate_row
        && mate_cells == 0
        && !pairhmm_keys.contains(&(MATE_QNAME.to_string(), MATE_FLAGS))
        && rust_by_key.len() == 122
        && java_only == 0
        && rust_only == 0
        && site.n_reads == 122;
    let matrix_ok = compared == 366 && differing == 0 && max_ulp == 0;
    let geometry_ok = (trim.padded_start, trim.padded_end) == JAVA_SPAN
        && (applied.trim_start, applied.trim_end) == JAVA_SPAN
        && read_repr == "51H97M @ 29455560 len=97"
        && ref_h.len == 169
        && ref_h.cigar == "169M"
        && lens == BTreeSet::from([169, 170, 172, 174])
        && hap.n == 24;
    let order_effect = gl_ulp == [0, 0, 0, 2, 0, 1]
        && replay.len() == 6
        && replay
            .iter()
            .zip(java_gl.iter())
            .all(|(r, j)| r.to_bits() == j.to_bits());
    let classification =
        if membership_ok && matrix_ok && geometry_ok && order_is_old && order_effect {
            "FAILED_MATE_PRODUCTION_PATCH_REPRODUCES_JAVA_MEMBERSHIP"
        } else if membership_ok {
            "FAILED_MATE_PRODUCTION_PATCH_PARTIAL"
        } else if rust_by_key.len() != 122 || java_only != 0 || rust_only != 0 {
            "FAILED_MATE_PRODUCTION_PATCH_REGRESSION"
        } else {
            "FAILED_MATE_PRODUCTION_PATCH_BLOCKED"
        };

    kv("active_region_reads", active_before.to_string());
    kv("ordinary_passing_reads", passing_filter.to_string());
    kv("failed_mate_exclusions", mate_fail_noted.to_string());
    kv("pairhmm_evidence", pairhmm_keys.len().to_string());
    kv("genotyping_evidence", site.n_reads.to_string());
    kv(
        "canonical_zero_row",
        if mate_row || mate_cells > 0 {
            "PRESENT"
        } else {
            "ABSENT"
        },
    );
    kv("java_only", java_only.to_string());
    kv("rust_only", rust_only.to_string());
    kv(
        "shared",
        rust_by_key.len().saturating_sub(rust_only).to_string(),
    );
    kv("compared_cells", compared.to_string());
    kv("differing_cells", differing.to_string());
    kv("max_abs", fmt_f(max_abs));
    kv("max_ulp", max_ulp.to_string());
    kv("final_gl", fmt_vec(&rust_gl));
    kv("final_pl", fmt_vec(&rust_pl));
    kv("gl_ulp", format!("{gl_ulp:?}"));
    kv("integer_pl", &integer);
    kv("emitted_subset_pl", &emitted_pl);
    kv(
        "rust_first",
        site.ad_row_qname.first().cloned().unwrap_or_default(),
    );
    kv("classification", classification);

    assert!(membership_ok, "{classification}");
    assert!(matrix_ok, "matrix differing={differing} ulp={max_ulp}");
    assert!(geometry_ok, "STR geometry");
    assert!(order_is_old, "read order changed");
    assert!(order_effect, "GL ULP {gl_ulp:?}");
    assert_eq!(
        classification,
        "FAILED_MATE_PRODUCTION_PATCH_REPRODUCES_JAVA_MEMBERSHIP"
    );
}
