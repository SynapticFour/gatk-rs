//! 6R.196: QUAL at `2:92307359 CT/C` uses Java length-based AF alt prior.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! Same PL-roundtrip GLs as Java; first split was SNP vs indel Dirichlet weight.
//! QD still QUAL/getDepth (6R.195). RankSums stay closed.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r196_qual_indel_af_prior -- --nocapture --test-threads=1
//! HOLDOUT_6R196=1 cargo test -p gatk-haplotypecaller --test holdout_6r196_qual_indel_af_prior -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const GAP_INTERVAL: &str = "2:92305500-92305850";
const TG_INTERVAL: &str = "2:92307200-92307550";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_307_359;
const CLOSED_GT: u64 = 92_305_634;
const CLOSED_AG: u64 = 92_305_635;
const CLOSED_AC: u64 = 92_305_716;
const CLOSED_INDEL: u64 = 92_307_324;
const CLOSED_TG: u64 = 92_307_333;
const CLOSED_MID: u64 = 92_316_347;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R196\t{key}\t{}", value.as_ref());
}

fn unique_indices(likelihoods: &[gatk_haplotypecaller::RegionReadLikelihood]) -> BTreeSet<usize> {
    likelihoods.iter().map(|c| c.read_index.get()).collect()
}

fn info_i32(info: &[InfoValue], key: &str) -> Option<i32> {
    for v in info {
        if let InfoValue::Integer(k, xs) = v {
            if k == key {
                return xs.first().copied();
            }
        }
    }
    None
}

fn info_f64(info: &[InfoValue], key: &str) -> Option<f64> {
    for v in info {
        if let InfoValue::Float(k, xs) = v {
            if k == key {
                return xs.first().copied();
            }
        }
    }
    None
}

fn info_has(info: &[InfoValue], key: &str) -> bool {
    info.iter().any(|v| match v {
        InfoValue::Flag(k)
        | InfoValue::Integer(k, _)
        | InfoValue::Float(k, _)
        | InfoValue::String(k, _)
        | InfoValue::Character(k, _) => k == key,
    })
}

#[test]
fn forensic_6r196_source_contracts() {
    kv("java_pin", JAVA_PIN);
    kv(
        "classification",
        "D — WRONG PROBABILITY / QUAL FORMULA (biallelic AF alt prior)",
    );
    kv(
        "production_change",
        "length-based AF alt Dirichlet weight at QUAL",
    );
    kv("target", "2:92307359 CT/C");

    let ann = include_str!("../src/variant_site_hc_annotations.rs");
    let af = include_str!("../src/af_calc.rs");
    let emit_gates = include_str!("../src/emit_gates.rs");
    let pipe = include_str!("../src/hc_genotyping_engine/genotype_site_pipeline.rs");
    let early = include_str!("../src/hc_genotyping_engine/genotype_site_early_template.rs");

    let annotate = ann
        .split("pub fn annotate_hc_variant_site")
        .nth(1)
        .expect("annotate");
    let annotate_body = annotate
        .split("fn read_offset_evidence_at_site")
        .next()
        .expect("body");
    assert!(
        annotate_body.contains("biallelic_alt_pseudocount(ref_allele, alt_allele")
            && annotate_body.contains("calculate_biallelic_af_em_with_alt_pseudocount"),
        "VCF QUAL AF must use Java length-based alt prior, not always snp_pseudocount"
    );
    assert!(
        !ann.contains("92307359") && !af.contains("92307359"),
        "no locus-specific QUAL patch"
    );
    assert!(
        af.contains("pub fn biallelic_alt_pseudocount")
            && af.contains("alt_allele.len() == ref_allele.len()")
            && af.contains("config.indel_pseudocount"),
        "Java: alt.length()==refLength → snpPseudocount else indelPseudocount"
    );
    let n2 = af
        .split("if n_alleles == 2 && span_del.is_none()")
        .nth(1)
        .expect("n2");
    assert!(
        n2.contains("biallelic_alt_pseudocount(alleles[0], alleles[1], config)"),
        "n==2 AF shortcut must not force the SNP prior"
    );
    assert!(
        emit_gates.contains("calculate_biallelic_af_em(")
            && !emit_gates.contains("biallelic_alt_pseudocount"),
        "emit threshold stays on default SNP-prior biallelic AF (emission predicates unchanged)"
    );
    assert!(
        pipe.contains("let annotation = annotation_likelihoods_from_stored_haplotypes")
            && early.contains("if is_cluster_tg_snp(&event)"),
        "6R.189 / 6R.186 paths stay"
    );
}

#[test]
fn forensic_6r196_qual_indel_af_prior_live() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
        return;
    }
    kv("java_pin", JAVA_PIN);
    kv("target", "2:92307359 CT/C");
    kv(
        "java_raw_qual",
        "31.60158353279359 (indel prior on PL 39,0,39; not inferred from printed 31.60)",
    );
    kv("rust_pre_fix_qual", "31.639400251058753 (SNP prior)");

    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let gap_specs = parse_intervals_cli_string(&dict, GAP_INTERVAL).expect("gap");
    let tg_specs = parse_intervals_cli_string(&dict, TG_INTERVAL).expect("tg");
    let gap_walk = traverse_assembly_region_walker(
        &dict,
        &gap_specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("gap walk");
    let tg_walk = traverse_assembly_region_walker(
        &dict,
        &tg_specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("tg walk");
    let gap_regions = flatten_assembly_regions(&gap_walk);
    let tg_regions = flatten_assembly_regions(&tg_walk);
    let covering_ag = gap_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_AG
                && r.end.get() >= CLOSED_AG
        })
        .expect("covering A/G");
    let covering_ac = gap_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_AC
                && r.end.get() >= CLOSED_AC
        })
        .expect("covering A/C");
    let covering_tg = tg_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("tg covering");

    let ag_outcome = HaplotypeCallerEngine::call_region(
        covering_ag,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("ag call")
    .expect("ag outcome");
    let ac_outcome = HaplotypeCallerEngine::call_region(
        covering_ac,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("ac call")
    .expect("ac outcome");
    let tg_outcome = HaplotypeCallerEngine::call_region(
        covering_tg,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("tg call")
    .expect("tg outcome");

    let closed_ag = ag_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_AG)
                && c.event.ref_allele == "A"
                && c.event.alt_allele == "G"
        })
        .expect("A/G");
    let closed_ac = ac_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_AC)
                && c.event.ref_allele == "A"
                && c.event.alt_allele == "C"
        })
        .expect("A/C");
    let target = tg_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(TARGET)
                && c.event.ref_allele == "CT"
                && c.event.alt_allele == "C"
        })
        .expect("CT/C");

    assert_eq!(unique_indices(&closed_ag.annotation_likelihoods).len(), 3);
    assert_eq!(unique_indices(&closed_ac.annotation_likelihoods).len(), 4);
    assert!(
        target.extra_alt_alleles.is_empty(),
        "biallelic emitted indel; QUAL is not a copied merged SPAN_DEL object"
    );
    assert!(
        target.qual_log10_p_error.is_none(),
        "QUAL comes from annotate AF on emitted GLs, not merged log10PError"
    );
    assert_eq!(target.genotype.format.pl_as_i32(), vec![39, 0, 39]);
    kv(
        "source_object",
        "PL-roundtrip GLs of the emitted biallelic genotype; AFCalculator P(no variant)",
    );
    kv("allele_set", "CT,C (REF len 2 ≠ ALT len 1 → indel prior)");

    let ag_emitted = try_emit_call_region_variants(
        covering_ag,
        &ag_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("ag emit");
    let ac_emitted = try_emit_call_region_variants(
        covering_ac,
        &ac_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("ac emit");
    let tg_emitted = try_emit_call_region_variants(
        covering_tg,
        &tg_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("tg emit");
    let gt_rec = ag_emitted
        .iter()
        .find(|r| r.position == CLOSED_GT && r.reference == "G")
        .expect("G/T");
    let ag_rec = ag_emitted
        .iter()
        .find(|r| r.position == CLOSED_AG && r.reference == "A")
        .expect("A/G");
    let ac_rec = ac_emitted
        .iter()
        .find(|r| r.position == CLOSED_AC && r.reference == "A")
        .expect("A/C");
    let indel_rec = tg_emitted
        .iter()
        .find(|r| r.position == CLOSED_INDEL && r.reference == "TTC")
        .expect("TTC/T");
    let tg_rec = tg_emitted
        .iter()
        .find(|r| r.position == CLOSED_TG && r.reference == "T")
        .expect("T/G");
    let target_rec = tg_emitted
        .iter()
        .find(|r| r.position == TARGET && r.reference == "CT")
        .expect("CT/C");

    assert_eq!(info_i32(&gt_rec.info, "DP"), Some(3));
    assert_eq!(info_i32(&ag_rec.info, "DP"), Some(3));
    assert_eq!(format!("{:.2}", ag_rec.quality.unwrap_or(0.0)), "78.32");
    assert!(!info_has(&ac_rec.info, "ReadPosRankSum"));
    assert_eq!(info_i32(&indel_rec.info, "DP"), Some(1));
    assert_eq!(format!("{:.2}", indel_rec.quality.unwrap_or(0.0)), "35.44");
    assert_eq!(info_i32(&tg_rec.info, "DP"), Some(1));
    assert_eq!(format!("{:.2}", tg_rec.quality.unwrap_or(0.0)), "35.48");

    let sample = target_rec.samples.first().expect("sample");
    assert_eq!(
        sample.gt.as_ref().map(|g| g.to_string()).as_deref(),
        Some("0/1")
    );
    assert_eq!(sample.ad.as_deref(), Some(&[1u32, 1][..]));
    assert_eq!(sample.dp.map(|d| d as i32), Some(2));
    assert_eq!(sample.gq.map(|g| g as i32), Some(39));
    let qual = target_rec.quality.expect("QUAL");
    kv("rust_qual_unrounded", format!("{qual}"));
    kv("rust_emitted_qual", format!("{qual:.2}"));
    assert!((qual - 31.60158353279359).abs() < 1e-9);
    assert_eq!(format!("{qual:.2}"), "31.60");
    let qd = info_f64(&target_rec.info, "QD").expect("QD");
    assert!((qd - qual / 2.0).abs() < 1e-12);
    assert_eq!(format!("{qd:.2}"), "15.80");
    assert!(info_has(&target_rec.info, "BaseQRankSum"));
    assert!(info_has(&target_rec.info, "MQRankSum"));
    assert!(info_has(&target_rec.info, "ReadPosRankSum"));

    let mid_specs = parse_intervals_cli_string(&dict, "2:92316200-92316580").expect("mid");
    let mid_walk = traverse_assembly_region_walker(
        &dict,
        &mid_specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("mid walk");
    let mid_regions = flatten_assembly_regions(&mid_walk);
    let covering_mid = mid_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_MID
                && r.end.get() >= CLOSED_MID
        })
        .expect("covering mid");
    let mid_outcome = HaplotypeCallerEngine::call_region(
        covering_mid,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("mid call")
    .expect("mid outcome");
    let mid_emitted = try_emit_call_region_variants(
        covering_mid,
        &mid_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("mid emit");
    let mid_rec = mid_emitted
        .iter()
        .find(|r| r.position == CLOSED_MID && r.reference == "G")
        .expect("G/A");
    assert_eq!(info_i32(&mid_rec.info, "DP"), Some(3));
    assert!((info_f64(&mid_rec.info, "MQ").unwrap_or(-1.0) - 40.25).abs() < 0.005);
}
