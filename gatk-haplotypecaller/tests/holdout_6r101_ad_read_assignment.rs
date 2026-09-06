//! 6R.101 holdout: first read-level FORMAT/AD divergence at the canonical T/C site.
//!
//! Skipped unless `HOLDOUT_6R101=1`. Coordinate-free contract lives in
//! `forensic_6r101_ad_read_assignment_contract`.
//!
//! Java live `DepthPerAlleleBySample.annotateWithLikelihoods` (seq 52 of
//! `HcParityAdAnnotationDump`) is 60×4 `TG,*,T,CG`, remaining `TG,CG`.
//! Four evidence rows have a 4-way unused best (`*` or `T`) and a remaining
//! remarg remaining-allele vote. That is Java's extra 2 TG + 2 CG versus the
//! unused-ALT permutation of 4-way counts.
//!
//! ```text
//! HOLDOUT_6R101=1 cargo test -p gatk-haplotypecaller --test holdout_6r101_ad_read_assignment -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::read_realignment::LOG_10_INFORMATIVE_THRESHOLD;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, take_colocated_merge_numerics,
    traverse_assembly_region_walker, try_emit_call_region_variants, AssemblyRegionCallDisposition,
    CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use serde_json::json;
use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const POS_SNP: u64 = 29_456_344;

#[derive(Clone, Debug)]
#[allow(dead_code)]
struct JavaRow {
    qname: String,
    flags: u16,
    four_best: String,
    four_delta: f64,
    four_informative: bool,
    four_counted: Option<String>,
    remarg_best: String,
    remarg_delta: f64,
    remarg_informative: bool,
    remarg_counted: Option<String>,
}

#[derive(Clone, Debug)]
#[allow(dead_code)]
struct RustRow {
    qname: String,
    flags: u16,
    four_best: String,
    four_delta: f64,
    four_informative: bool,
    four_counted: Option<String>,
    remarg_best: String,
    remarg_delta: f64,
    remarg_informative: bool,
    remarg_counted: Option<String>,
}

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn parse_counted(s: &str) -> Option<String> {
    if s.is_empty() || s == "." {
        None
    } else {
        Some(s.to_string())
    }
}

fn load_java_rows() -> Vec<JavaRow> {
    let path =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/6r101_java_seq52_ad_assignment.tsv");
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path:?}: {e}"));
    let mut rows = Vec::new();
    for (i, line) in text.lines().enumerate() {
        if i == 0 || line.is_empty() {
            continue;
        }
        let c: Vec<&str> = line.split('\t').collect();
        assert!(c.len() >= 17, "java tsv cols {}", c.len());
        rows.push(JavaRow {
            qname: c[1].to_string(),
            flags: c[2].parse().unwrap(),
            four_best: c[3].to_string(),
            four_delta: c[7].parse().unwrap(),
            four_informative: c[8] == "true",
            four_counted: parse_counted(c[9]),
            remarg_best: c[10].to_string(),
            remarg_delta: c[14].parse().unwrap(),
            remarg_informative: c[15] == "true",
            remarg_counted: parse_counted(c[16]),
        });
    }
    rows
}

fn vote(lls: &[f64], names: &[String]) -> (String, String, f64, f64, bool, Option<String>) {
    let mut best_i = 0usize;
    let mut best = f64::NEG_INFINITY;
    let mut second = f64::NEG_INFINITY;
    let mut second_i = 0usize;
    for (i, &ll) in lls.iter().enumerate() {
        let ll = if ll.is_finite() { ll } else { -50.0 };
        if ll > best {
            second = best;
            second_i = best_i;
            best = ll;
            best_i = i;
        } else if ll > second {
            second = ll;
            second_i = i;
        }
    }
    let delta = if best == second { 0.0 } else { best - second };
    let inf = best.is_finite() && delta.abs() > LOG_10_INFORMATIVE_THRESHOLD;
    let best_n = names
        .get(best_i)
        .cloned()
        .unwrap_or_else(|| best_i.to_string());
    let second_n = names
        .get(second_i)
        .cloned()
        .unwrap_or_else(|| second_i.to_string());
    let counted = if inf { Some(best_n.clone()) } else { None };
    (best_n, second_n, best, delta, inf, counted)
}

fn remaining_names(long_ref: &str, alts: &[String], keep: &[usize]) -> Vec<String> {
    let all: Vec<String> = std::iter::once(long_ref.to_string())
        .chain(alts.iter().cloned())
        .collect();
    keep.iter().filter_map(|&i| all.get(i).cloned()).collect()
}

fn rust_rows_from_snap(snap: &gatk_haplotypecaller::ColocatedMergeNumerics) -> Vec<RustRow> {
    let names: Vec<String> = std::iter::once(snap.long_ref.clone())
        .chain(snap.alts.iter().cloned())
        .collect();
    let rem_names = remaining_names(&snap.long_ref, &snap.alts, &snap.remaining_keep_indices);
    snap.ad_row_qname
        .iter()
        .enumerate()
        .map(|(i, qn)| {
            let lls = &snap.ad_row_lls[i];
            let (four_best, _, _, four_delta, four_inf, four_counted) = vote(lls, &names);
            let rem_lls: Vec<f64> = snap
                .remaining_keep_indices
                .iter()
                .filter_map(|&k| lls.get(k).copied())
                .collect();
            let (remarg_best, _, _, remarg_delta, rem_inf, remarg_counted) =
                vote(&rem_lls, &rem_names);
            RustRow {
                qname: qn.clone(),
                flags: snap.ad_row_flags[i],
                four_best,
                four_delta,
                four_informative: four_inf,
                four_counted,
                remarg_best,
                remarg_delta,
                remarg_informative: rem_inf,
                remarg_counted,
            }
        })
        .collect()
}

fn counted(c: &Option<String>) -> &str {
    c.as_deref().unwrap_or("UNINFORMATIVE")
}

#[test]
fn holdout_6r101_first_ad_assignment_divergence() {
    if std::env::var("HOLDOUT_6R101").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R101=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    assert!(ref_fasta.is_file() && bam.is_file());

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
    let covering: Vec<_> = regions
        .iter()
        .filter(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= POS_SNP
                && r.end.get() >= POS_SNP
        })
        .collect();
    assert_eq!(covering.len(), 1);
    let outcome = HaplotypeCallerEngine::call_region(
        covering[0],
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let emitted = try_emit_call_region_variants(
        covering[0],
        &outcome,
        "SAMPLE",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .unwrap_or_default();
    let vcf = emitted
        .iter()
        .find(|r| {
            r.position == POS_SNP && r.reference == "T" && r.alternate.iter().any(|a| a == "C")
        })
        .expect("T/C");
    let snap = take_colocated_merge_numerics()
        .into_iter()
        .find(|n| n.loc == POS_SNP)
        .expect("snap");

    let java = load_java_rows();
    let rust = rust_rows_from_snap(&snap);
    assert_eq!(java.len(), 60);
    assert_eq!(rust.len(), snap.n_reads);
    assert_eq!(snap.n_reads, 60);

    let java_keys: HashSet<(String, u16)> =
        java.iter().map(|r| (r.qname.clone(), r.flags)).collect();
    let rust_keys: HashSet<(String, u16)> =
        rust.iter().map(|r| (r.qname.clone(), r.flags)).collect();
    let java_only: Vec<_> = java_keys.difference(&rust_keys).collect();
    let rust_only: Vec<_> = rust_keys.difference(&java_keys).collect();
    eprintln!(
        "6R.101 membership java={} rust={} JAVA_ONLY={} RUST_ONLY={}",
        java.len(),
        rust.len(),
        java_only.len(),
        rust_only.len()
    );
    assert!(
        java_only.is_empty() && rust_only.is_empty(),
        "AD_EVIDENCE_MEMBERSHIP: {java_only:?} {rust_only:?}"
    );

    let rust_by: HashMap<(String, u16), &RustRow> = rust
        .iter()
        .map(|r| ((r.qname.clone(), r.flags), r))
        .collect();

    let mut java_inf = 0usize;
    let mut rust_inf = 0usize;
    let mut java_tg = 0usize;
    let mut rust_tg = 0usize;
    let mut java_cg = 0usize;
    let mut rust_cg = 0usize;
    let mut java_uninf = 0usize;
    let mut rust_uninf = 0usize;
    let mut java_only_inf = 0usize;
    let mut rust_only_inf = 0usize;
    let mut java_only_tg = 0usize;
    let mut rust_only_tg = 0usize;
    let mut java_only_cg = 0usize;
    let mut rust_only_cg = 0usize;
    let mut diverge = Vec::new();

    for j in &java {
        let r = rust_by.get(&(j.qname.clone(), j.flags)).expect("matched");
        if j.remarg_informative {
            java_inf += 1;
        } else {
            java_uninf += 1;
        }
        if r.remarg_informative {
            rust_inf += 1;
        } else {
            rust_uninf += 1;
        }
        match counted(&j.remarg_counted) {
            "TG" => java_tg += 1,
            "CG" => java_cg += 1,
            _ => {}
        }
        match counted(&r.remarg_counted) {
            "TG" => rust_tg += 1,
            "CG" => rust_cg += 1,
            _ => {}
        }
        if j.remarg_informative && !r.remarg_informative {
            java_only_inf += 1;
        }
        if r.remarg_informative && !j.remarg_informative {
            rust_only_inf += 1;
        }
        if counted(&j.remarg_counted) == "TG" && counted(&r.remarg_counted) != "TG" {
            java_only_tg += 1;
        }
        if counted(&r.remarg_counted) == "TG" && counted(&j.remarg_counted) != "TG" {
            rust_only_tg += 1;
        }
        if counted(&j.remarg_counted) == "CG" && counted(&r.remarg_counted) != "CG" {
            java_only_cg += 1;
        }
        if counted(&r.remarg_counted) == "CG" && counted(&j.remarg_counted) != "CG" {
            rust_only_cg += 1;
        }

        let remarg_match = counted(&j.remarg_counted) == counted(&r.remarg_counted);
        let four_match = counted(&j.four_counted) == counted(&r.four_counted);
        if !remarg_match || counted(&j.four_counted) != counted(&j.remarg_counted) {
            let first = if !four_match && j.four_best != r.four_best {
                "AD_BEST_ALLELE_SELECTION"
            } else if j.remarg_informative != r.remarg_informative && j.remarg_best == r.remarg_best
            {
                "AD_INFORMATIVENESS_PREDICATE"
            } else if counted(&j.four_counted) != counted(&j.remarg_counted) {
                "AD_BEST_ALLELE_SELECTION"
            } else {
                "OTHER"
            };
            diverge.push(json!({
                "qname": j.qname,
                "flags": j.flags,
                "java_best": j.remarg_best,
                "rust_best": r.remarg_best,
                "java_delta": j.remarg_delta,
                "rust_delta": r.remarg_delta,
                "java_informative": j.remarg_informative,
                "rust_informative": r.remarg_informative,
                "counted_java": counted(&j.remarg_counted),
                "counted_rust_remarg": counted(&r.remarg_counted),
                "java_four_counted": counted(&j.four_counted),
                "rust_four_counted": counted(&r.four_counted),
                "java_four_best": j.four_best,
                "rust_four_best": r.four_best,
                "java_four_delta": j.four_delta,
                "first_divergence": first,
            }));
        }
        let _ = remarg_match;
    }

    let permute_div: Vec<_> = java
        .iter()
        .filter(|j| counted(&j.four_counted) != counted(&j.remarg_counted))
        .map(|j| {
            let r = rust_by[&(j.qname.clone(), j.flags)];
            (
                j.qname.clone(),
                j.flags,
                j.four_best.clone(),
                r.four_best.clone(),
                counted(&j.four_counted).to_string(),
                counted(&j.remarg_counted).to_string(),
                counted(&r.remarg_counted).to_string(),
            )
        })
        .collect();

    let vcf_ad = vcf.samples[0].ad.clone().unwrap_or_default();
    let vcf_pl = vcf.samples[0].pl.clone().unwrap_or_default();
    let extra_remaining_ad: Vec<_> = permute_div
        .iter()
        .filter(|(_, _, _, _, _, jrem, _)| jrem == "TG" || jrem == "CG")
        .cloned()
        .collect();
    let unused_then_uninf = permute_div.len() - extra_remaining_ad.len();

    let doc = json!({
        "java_annotation": {
            "evidence": 60,
            "alleles": ["TG", "*", "T", "CG"],
            "remaining": ["TG", "CG"],
            "informative": java_inf,
            "uninformative": java_uninf,
            "TG": java_tg,
            "CG": java_cg,
            "ad": [36, 19],
        },
        "rust_annotation": {
            "evidence": rust.len(),
            "informative": rust_inf,
            "uninformative": rust_uninf,
            "TG": rust_tg,
            "CG": rust_cg,
            "remarg": snap.subset_ad_remarginalized,
            "permute": snap.subset_ad_permuted,
            "vcf_ad": vcf_ad,
        },
        "JAVA_ONLY_INFORMATIVE": java_only_inf,
        "RUST_ONLY_INFORMATIVE": rust_only_inf,
        "JAVA_ONLY_TG": java_only_tg,
        "RUST_ONLY_TG": rust_only_tg,
        "JAVA_ONLY_CG": java_only_cg,
        "RUST_ONLY_CG": rust_only_cg,
        "four_way_vs_remarg_reads": permute_div,
        "extra_remaining_ad_reads": extra_remaining_ad,
        "unused_four_way_then_uninformative_remarg": unused_then_uninf,
        "divergent_remarg_assignment": diverge,
        "first_divergent_operation": "AD_BEST_ALLELE_SELECTION",
        "note": "same 60 reads; remarg remaining matches; FORMAT AD was 4-way permute of unused ALTs",
        "vcf": {
            "gt": vcf.samples[0].gt.as_ref().map(|g| g.alleles.clone()),
            "ad": vcf_ad,
            "pl": vcf_pl,
            "qual": vcf.quality,
        },
        "java_oracle": {"gt": [0, 1], "ad": [36, 19], "pl": [542, 0, 1353], "qual": 510.06},
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());

    assert_eq!(java_only_inf, 0);
    assert_eq!(rust_only_inf, 0);
    assert_eq!(java_only_tg, 0);
    assert_eq!(rust_only_tg, 0);
    assert_eq!(java_only_cg, 0);
    assert_eq!(rust_only_cg, 0);
    assert_eq!(java_tg, 36);
    assert_eq!(java_cg, 19);
    assert_eq!(rust_tg, 36);
    assert_eq!(rust_cg, 19);
    assert_eq!(snap.subset_ad_remarginalized, vec![36, 19]);
    assert_eq!(snap.subset_ad_permuted, vec![34, 17]);
    assert_eq!(
        extra_remaining_ad.len(),
        4,
        "Java extra 2 TG + 2 CG vs 4-way permute: {extra_remaining_ad:?}"
    );
    assert_eq!(extra_remaining_ad.iter().filter(|t| t.5 == "TG").count(), 2);
    assert_eq!(extra_remaining_ad.iter().filter(|t| t.5 == "CG").count(), 2);
    for (qn, flags, j4, r4, four_c, jrem, rrem) in &permute_div {
        eprintln!(
            "6R.101 FOUR_VS_REMARG {qn} flags={flags} java_four={j4} rust_four={r4} four_counted={four_c} java_remarg={jrem} rust_remarg={rrem}"
        );
        assert_eq!(j4, r4, "4-way best allele already matches for {qn}");
        assert_eq!(jrem, rrem, "remaining remarg already matches for {qn}");
        assert_ne!(four_c.as_str(), jrem.as_str());
        assert!(j4 == "*" || j4 == "T");
    }
    assert_eq!(vcf_pl, vec![542u32, 0, 1353]);
    assert_eq!(
        vcf_ad,
        vec![36u32, 19],
        "FORMAT AD is remaining remarg, not 4-way permute"
    );
    assert_eq!(
        vcf_ad.iter().map(|&x| x as i32).collect::<Vec<_>>(),
        snap.subset_ad_remarginalized
    );
    assert_ne!(
        snap.subset_ad_permuted, snap.subset_ad_remarginalized,
        "permute remains a distinct unused-ALT slice"
    );
}
