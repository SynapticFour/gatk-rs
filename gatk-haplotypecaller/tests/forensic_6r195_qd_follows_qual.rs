//! 6R.195: proof-only. QD at `2:92307359 CT/C` follows QUAL / AD-depth.
//! PRODUCTION CHANGE: NONE.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! RankSums stay closed. QUAL is an upstream frozen residual.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r195_qd_follows_qual -- --nocapture --test-threads=1
//! HOLDOUT_6R195=1 cargo test -p gatk-haplotypecaller --test holdout_6r195_qd_follows_qual -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::{Genotype, InfoValue};
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::genotyping::GenotypeFormatFields;
use gatk_haplotypecaller::variant_site_hc_annotations::qd_depth_for_variant;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const GAP_INTERVAL: &str = "2:92305500-92305850";
const TG_INTERVAL: &str = "2:92307200-92307550";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_VCF_REL: &str = "parity/reports/6r43/p12_indel_mix/java.vcf";
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
    eprintln!("6R195\t{key}\t{}", value.as_ref());
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

fn parse_java_target(root: &Path) -> Option<(String, String)> {
    let text = fs::read_to_string(root.join(JAVA_VCF_REL)).ok()?;
    for line in text.lines() {
        if line.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 8 {
            continue;
        }
        if cols[0] != "2" || cols[1] != "92307359" || cols[3] != "CT" {
            continue;
        }
        let qd = cols[7]
            .split(';')
            .find_map(|kv| kv.strip_prefix("QD=").map(str::to_string))?;
        return Some((cols[5].to_string(), qd));
    }
    None
}

#[test]
fn forensic_6r195_source_contracts() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv(
        "classification",
        "E — NO QD DIVERGENCE: QD correctly follows an already-divergent QUAL",
    );
    kv("target", "2:92307359 CT/C");

    let ann = include_str!("../src/variant_site_hc_annotations.rs");
    let plugin = include_str!("../src/annotator/plugins/qual_by_depth.rs");
    let emit = include_str!("../src/region_vcf_emit.rs");
    let vcf = include_str!("../../gatk-core/src/io/vcf.rs");
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
        annotate_body.contains("let qd_depth = qd_depth_for_variant")
            && annotate_body.contains("qual_by_depth::raw_qual_by_depth(qual, qd_depth)"),
        "QD numerator is the same unrounded QUAL used for the QUAL column; denominator is getDepth"
    );
    assert!(
        !annotate_body.contains("coverage_evidence_count")
            || annotate_body.contains("6R.174: INFO DP is Java `Coverage.evidenceCount`"),
        "INFO DP coverage remains a separate field from QD depth"
    );
    let qd_depth_fn = ann
        .split("pub fn qd_depth_for_variant")
        .nth(1)
        .expect("qd_depth")
        .split("fn is_het_or_hom_var")
        .next()
        .expect("body");
    assert!(
        qd_depth_fn.contains("fields.ad") && qd_depth_fn.contains("alt_ad > 1"),
        "Java QualByDepth.getDepth: AD sum for het/hom-var; ADrestricted only if alt AD > 1"
    );
    assert!(
        !qd_depth_fn.contains("92307359") && !ann.contains("92307359"),
        "no locus-specific QD patch"
    );

    assert!(
        plugin.contains("pub fn raw_qual_by_depth(qual: f64, dp: i32) -> f64")
            && plugin.contains("qual / (dp as f64)")
            && plugin.contains("MAX_QD_BEFORE_FIXING: f64 = 35.0")
            && plugin.contains("pub fn apply_fix_too_high_qd_to_vcf_records"),
        "raw QD is QUAL/depth; jitter only if QD >= 35, applied later in emit order"
    );
    assert!(
        emit.contains("quality: if hom_ref { None } else { Some(ann.qual) }")
            && emit.contains("InfoValue::Float(\"QD\".to_string(), vec![ann.qd])"),
        "emitted QUAL and QD share the annotation bundle; QD is not recomputed from the printed QUAL string"
    );
    assert!(
        vcf.contains("|q| format!(\"{:.2}\", q)") || vcf.contains("format!(\"{:.2}\", q)"),
        "QUAL column is HTSJDK-style %.2f of the unrounded phred"
    );
    assert!(
        pipe.contains("let annotation = annotation_likelihoods_from_stored_haplotypes")
            && early.contains("if is_cluster_tg_snp(&event)"),
        "6R.189 / 6R.186 paths stay"
    );
}

#[test]
fn forensic_6r195_java_getdepth_ad_one_one() {
    kv("java_path", "QualByDepth.annotate → -10*vc.getLog10PError() / getDepth(genotypes, likelihoods) → fixTooHighQD if QD>=35 → String.format(%.2f)");
    let gt = Genotype {
        alleles: vec![0, 1],
        phased: false,
    };
    let fields = GenotypeFormatFields::from_wire(vec![39, 0, 39], 39, vec![1, 1], 2);
    let depth = qd_depth_for_variant(&gt, &fields);
    kv("java_getdepth_ad_1_1", depth.to_string());
    assert_eq!(
        depth, 2,
        "AD=1,1 het: total AD=2, alt AD=1 is not >1 so ADrestricted unused; depth=2"
    );
    let dummy_info_dp = 99;
    assert_ne!(
        depth, dummy_info_dp,
        "QD depth is FORMAT AD sum, not an arbitrary INFO DP"
    );
}

#[test]
fn forensic_6r195_qd_follows_qual() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
        return;
    }
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv("target", "2:92307359 CT/C");

    if let Some((j_qual, j_qd)) = parse_java_target(&root) {
        kv("java_emitted_qual", &j_qual);
        kv("java_emitted_qd", &j_qd);
        assert_eq!(j_qual, "31.60");
        assert_eq!(j_qd, "15.80");
        let jq: f64 = j_qual.parse().expect("java qual");
        assert_eq!(
            format!("{:.2}", jq / 2.0),
            j_qd,
            "Java printed QD is printed QUAL / 2 at %.2f (same displayed QUAL bits)"
        );
    } else {
        kv(
            "java_vcf",
            "missing p12_indel_mix java.vcf; live Rust proof still runs",
        );
    }

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
    let closed_tg = tg_outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_TG)
                && c.event.ref_allele == "T"
                && c.event.alt_allele == "G"
        })
        .expect("T/G");
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
    assert_eq!(unique_indices(&closed_tg.annotation_likelihoods).len(), 1);
    assert_eq!(target.genotype.format.ad_as_i32(), vec![1, 1]);
    assert_eq!(target.genotype.format.dp.as_i32(), 2);

    let gt = Genotype {
        alleles: vec![0, 1],
        phased: false,
    };
    let qd_depth = qd_depth_for_variant(&gt, &target.genotype.format);
    kv("rust_qd_depth", qd_depth.to_string());
    assert_eq!(qd_depth, 2);

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
    assert!((info_f64(&ag_rec.info, "MQ").unwrap_or(-1.0) - 41.96).abs() < 0.005);
    assert!(!info_has(&ac_rec.info, "ReadPosRankSum"));
    assert_eq!(info_i32(&indel_rec.info, "DP"), Some(1));
    assert_eq!(info_i32(&tg_rec.info, "DP"), Some(1));

    let sample = target_rec.samples.first().expect("sample");
    assert_eq!(
        sample.gt.as_ref().map(|g| g.to_string()).as_deref(),
        Some("0/1")
    );
    assert_eq!(sample.ad.as_deref(), Some(&[1u32, 1][..]));
    assert_eq!(sample.dp.map(|d| d as i32), Some(2));
    assert_eq!(sample.gq.map(|g| g as i32), Some(39));

    let qual = target_rec.quality.expect("QUAL");
    let qd = info_f64(&target_rec.info, "QD").expect("QD");
    kv("rust_qual_unrounded", format!("{qual}"));
    kv("rust_emitted_qual_2dp", format!("{qual:.2}"));
    kv("rust_qd_raw", format!("{qd}"));
    kv("rust_qd_2dp", format!("{qd:.2}"));
    kv(
        "rust_qual_over_depth",
        format!("{}", qual / qd_depth as f64),
    );
    assert!((qual - 31.60).abs() < 0.05);
    assert_eq!(format!("{qual:.2}"), "31.60");
    assert!(
        (qd - qual / qd_depth as f64).abs() < 1e-12,
        "emitted QD must be the unrounded QUAL / getDepth, got qd={qd} qual/2={}",
        qual / 2.0
    );
    assert!(qd < 35.0, "fixTooHighQD is a no-op below 35");
    assert_eq!(format!("{qd:.2}"), "15.80");
    assert_eq!(
        format!("{:.2}", 31.60_f64 / 2.0),
        "15.80",
        "QD remains QUAL/2 (6R.195 E); 6R.196 indel AF prior closed the QUAL split"
    );
    assert!(info_has(&target_rec.info, "BaseQRankSum"));
    assert!(info_has(&target_rec.info, "MQRankSum"));
    assert!(info_has(&target_rec.info, "ReadPosRankSum"));
    kv(
        "classification",
        "E — NO QD DIVERGENCE: QD correctly follows an already-divergent QUAL",
    );

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
