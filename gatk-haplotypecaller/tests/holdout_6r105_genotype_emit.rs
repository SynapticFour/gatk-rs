//! 6R.105 live dump: genotype assignment / emit threshold at HOLDOUT_6R53
//! `20:29455388 C/T`. Skipped unless `HOLDOUT_6R105=1`.
//!
//! Original 6R.105 finding: given alleles `[C, T]`, Java `PL=0,6,1780` assigns
//! GT `0/0` and `calculateGenotypes` is null; the then-live Rust `PL=81,0,36`
//! assigned GT `0/1` and emitted. Emit predicates agree given the same GLs.
//! 6R.106/6R.107 proved the het PL was L9 SparsePlShape, not the calculator.
//! After HomAltStrong gating, live production matches Java: no genotyped C/T,
//! no emit. The GL→emit contract remains `PL=0,6,1780` → GT `0/0` → emit false.
//!
//! ```text
//! HOLDOUT_6R105=1 cargo test -p gatk-haplotypecaller --test holdout_6r105_genotype_emit -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::emit_gates::{
    java_emit_would_pass, passes_hc_variant_emit_biallelic, passes_java_emit_not_hom_ref,
};
use gatk_haplotypecaller::event_map::VariationEvent;
use gatk_haplotypecaller::genotyping::{
    best_pl_index, diploid_genotype_alleles_from_pl_index, emit_genotype_format_fields,
};
use gatk_haplotypecaller::hc_genotyping_engine::{
    java_emit_af_decision, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use serde_json::json;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/java.vcf";
const RUST_VCF_REL: &str = "parity/reports/6r43/chr20_tiny/rust.vcf";
const COVERING: (u64, u64) = (29_455_300, 29_455_559);
const TARGET: u64 = 29_455_388;
const JAVA_STAND_CALL_CONF: f64 = 30.0;

fn keep_alleles_from_assigned_gt(n_alleles: usize, gt: &[i32]) -> Vec<usize> {
    let mut keep = vec![0];
    for i in 1..n_alleles {
        if gt.iter().any(|&g| g >= 0 && (g as usize) == i) {
            keep.push(i);
        }
    }
    keep
}

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn vcf_has_pos(path: &Path, pos: u64) -> bool {
    for line in std::fs::read_to_string(path).unwrap_or_default().lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<_> = line.split('\t').collect();
        if f.len() < 5 || f[0] != "20" {
            continue;
        }
        if f[1].parse::<u64>().ok() == Some(pos) {
            return true;
        }
    }
    false
}

fn vcf_record(path: &Path, pos: u64) -> Option<serde_json::Value> {
    for line in std::fs::read_to_string(path).unwrap_or_default().lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<_> = line.split('\t').collect();
        if f.len() < 10 || f[0] != "20" {
            continue;
        }
        if f[1].parse::<u64>().ok()? != pos {
            continue;
        }
        let fmt: Vec<_> = f[8].split(':').collect();
        let samp: Vec<_> = f[9].split(':').collect();
        let get = |k: &str| {
            fmt.iter()
                .position(|x| *x == k)
                .and_then(|i| samp.get(i))
                .unwrap_or(&".")
                .to_string()
        };
        return Some(json!({
            "alleles": format!("{}/{}", f[3], f[4]),
            "qual": f[5],
            "gt": get("GT"),
            "ad": get("AD"),
            "pl": get("PL"),
            "gq": get("GQ"),
        }));
    }
    None
}

#[test]
fn holdout_6r105_genotype_emit_dump() {
    if std::env::var("HOLDOUT_6R105").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R105=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    let java_vcf = root.join(JAVA_VCF_REL);
    let rust_vcf = root.join(RUST_VCF_REL);
    assert!(ref_fasta.is_file(), "missing {}", ref_fasta.display());
    assert!(bam.is_file(), "missing {}", bam.display());

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
            ) && r.start.get() == COVERING.0
                && r.end.get() == COVERING.1
        })
        .expect("covering ActiveFull");
    let args = CallRegionArgs::strict_java();
    let outcome = HaplotypeCallerEngine::call_region(covering, &dict, &ref_fasta, &args)
        .expect("call")
        .expect("ActiveFull outcome");

    let calls: Vec<_> = outcome
        .genotyped_calls
        .iter()
        .filter(|c| c.event.start_1based.get() == TARGET)
        .collect();
    let event = VariationEvent::from_alleles("20", TARGET, "C", "T");
    let calculator_pl = [0, 6, 1780];
    let gls: Vec<f64> = calculator_pl.iter().map(|&p| (p as f64) / -10.0).collect();
    let format = emit_genotype_format_fields(&gls, &[44, 4]).expect("fmt");
    let pl = format.pl_as_i32();
    let gt = diploid_genotype_alleles_from_pl_index(2, best_pl_index(&format.pl));
    let keep = keep_alleles_from_assigned_gt(2, &gt);

    let af10 = java_emit_af_decision(&gls, DEFAULT_STAND_EMIT_CONFIDENCE).expect("af10");
    let af30 = java_emit_af_decision(&gls, JAVA_STAND_CALL_CONF).expect("af30");
    let not_hom_ref = passes_java_emit_not_hom_ref(&gls, &format);
    let emit_would_10 =
        java_emit_would_pass(&event, &gls, &format, DEFAULT_STAND_EMIT_CONFIDENCE, &[])
            .expect("emit10");
    let emit_would_30 =
        java_emit_would_pass(&event, &gls, &format, JAVA_STAND_CALL_CONF, &[]).expect("emit30");
    let site_af_10 =
        passes_hc_variant_emit_biallelic(&gls, DEFAULT_STAND_EMIT_CONFIDENCE).expect("site10");
    let site_af_30 = passes_hc_variant_emit_biallelic(&gls, JAVA_STAND_CALL_CONF).expect("site30");

    let emitted =
        try_emit_call_region_variants(covering, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let rust_emitted = emitted.iter().any(|r| r.position == TARGET);

    let doc = json!({
        "holdout": "20:29455388 C/T",
        "entering_alleles": format!("{}/{}", event.ref_allele, event.alt_allele),
        "live_genotyped": !calls.is_empty(),
        "n_gl": gls.len(),
        "log10_gl": gls,
        "pl": pl,
        "assigned_gt_from_pl": gt,
        "output_allele_keep_indices": keep,
        "stand_emit_confidence_rust": DEFAULT_STAND_EMIT_CONFIDENCE,
        "java_stand_call_conf": JAVA_STAND_CALL_CONF,
        "passes_java_emit_not_hom_ref": not_hom_ref,
        "passes_hc_variant_emit_biallelic_10": site_af_10,
        "passes_hc_variant_emit_biallelic_30": site_af_30,
        "java_emit_would_pass_10": emit_would_10,
        "java_emit_would_pass_30": emit_would_30,
        "af_decision_10": {
            "phred_scaled": af10.phred_scaled,
            "site_is_monomorphic": af10.site_is_monomorphic,
            "alt_plausible": af10.alt_plausible,
            "log10_posterior_no_variant": af10.log10_posterior_no_variant,
            "passes_emit": af10.passes_emit,
        },
        "af_decision_30": {
            "phred_scaled": af30.phred_scaled,
            "site_is_monomorphic": af30.site_is_monomorphic,
            "alt_plausible": af30.alt_plausible,
            "log10_posterior_no_variant": af30.log10_posterior_no_variant,
            "passes_emit": af30.passes_emit,
        },
        "rust_emitted": rust_emitted,
        "frozen_vcf": {
            "java_has_record": vcf_has_pos(&java_vcf, TARGET),
            "rust": vcf_record(&rust_vcf, TARGET),
        },
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());

    assert!(
        !vcf_has_pos(&java_vcf, TARGET),
        "Java VCF still has no record"
    );
    assert!(
        vcf_record(&rust_vcf, TARGET).is_none(),
        "6R.107: frozen Rust VCF must not contain C/T"
    );
    assert!(
        calls.is_empty(),
        "6R.107: Rust must not genotype C/T at 20:29455388"
    );
    assert!(!rust_emitted, "6R.107: Rust must not emit C/T");
    assert_eq!(pl.as_slice(), [0, 6, 1780]);
    assert_eq!(gt, vec![0, 0], "calculator PL assigns GT 0/0");
    assert_eq!(keep, vec![0], "hom-ref output-allele subset drops T");
    assert!(!not_hom_ref);
    assert!(!emit_would_10);
    assert!(!emit_would_30);
    assert!(af30.site_is_monomorphic);
    assert!(!af30.passes_emit);
    assert!(!site_af_10);
    assert!(!site_af_30);
}
