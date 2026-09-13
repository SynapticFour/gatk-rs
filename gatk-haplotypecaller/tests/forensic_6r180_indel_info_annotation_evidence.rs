//! 6R.180: INFO DP/MQ/SOR consume per-variant genotyping AlleleLikelihoods.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! Target `2:92307324 TTC/T`. Formulas unchanged. FORMAT/QUAL unchanged.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r180_indel_info_annotation_evidence -- --nocapture --test-threads=1
//! HOLDOUT_6R180=1 cargo test -p gatk-haplotypecaller --test holdout_6r180_indel_info -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::hc_genotyping_engine::java_alignment_read_overlaps_interval;
use gatk_haplotypecaller::variant_site_hc_annotations::{
    coverage_evidence_count, rms_mapping_quality_raw, rms_mapping_quality_sample_mapqs,
    strand_bias_contingency_table, HcStrandBiasLikelihoods, STRAND_ODDS_RATIO_MIN_COUNT,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92307200-92307550";
const CLOSED_SNP: u64 = 92_305_634;
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_307_324;
const MERGED_REF: &str = "TTC";
const MERGED_ALT: &str = "T";
const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const JAVA_READ: &str = "H06JUADXX130110:1:1101:10052:88682";
const FLAG_REVERSE: u16 = 0x10;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R180\t{key}\t{}", value.as_ref());
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

fn java_calculate_sor(ref_fw: u32, ref_rv: u32, alt_fw: u32, alt_rv: u32) -> f64 {
    let t00 = f64::from(ref_fw) + 1.0;
    let t01 = f64::from(ref_rv) + 1.0;
    let t10 = f64::from(alt_fw) + 1.0;
    let t11 = f64::from(alt_rv) + 1.0;
    let ratio = (t00 / t01) * (t11 / t10) + (t01 / t00) * (t10 / t11);
    let ref_ratio = t00.min(t01) / t00.max(t01);
    let alt_ratio = t10.min(t11) / t10.max(t11);
    ratio.ln() + ref_ratio.ln() - alt_ratio.ln()
}

#[test]
fn forensic_6r180_formulas_untouched() {
    kv("java_pin", JAVA_PIN);
    kv(
        "production_change",
        "per-variant annotation AlleleLikelihoods",
    );
    let ann = include_str!("../src/variant_site_hc_annotations.rs");
    let emit = include_str!("../src/region_vcf_emit.rs");
    let pipe = include_str!("../src/hc_genotyping_engine/genotype_site_pipeline.rs");
    let ge = include_str!("../src/hc_genotyping_engine/mod.rs");
    assert!(
        ann.contains("6R.174: INFO DP is Java `Coverage.evidenceCount`")
            && ann.contains("coverage_evidence_count(")
            && ann.contains("mq_rms_of_sample_evidence")
            && ann.contains("strand_bias_sample_counts"),
        "DP/MQ/SOR formulas must stay closed"
    );
    assert!(
        emit.contains("annotation_likelihoods_for_call")
            && emit.contains("6R.180: bind INFO DP/MQ/SOR"),
        "6R.180 shared evidence binding"
    );
    assert!(
        pipe.contains("with_annotation_likelihoods")
            && ge.contains("fn per_variant_annotation_likelihoods"),
        "genotyping subset is passed through, not recomputed at emit"
    );
    assert!(
        !ann.contains("92307324") && !emit.contains("92307324") && !pipe.contains("92307324"),
        "no locus pin"
    );
    let one = java_calculate_sor(0, 0, 0, 1);
    assert!((one - 1.6094379124341003).abs() < 1e-12);
    let mq = rms_mapping_quality_raw(&[44]);
    assert_eq!(mq.map(|t| t.2), Some(44.0));
}

#[test]
fn forensic_6r180_per_variant_annotation_evidence() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
        return;
    }
    kv("java_pin", JAVA_PIN);
    kv("target", "2:92307324 TTC/T");

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
        .expect("covering");
    let outcome = HaplotypeCallerEngine::call_region(
        covering,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("outcome");
    let call = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(TARGET)
                && c.event.ref_allele == MERGED_REF
                && c.event.alt_allele == MERGED_ALT
        })
        .expect("call");
    assert_eq!(call.genotype.format.dp.as_i32(), 1);
    assert_eq!(call.genotype.format.ad_as_i32(), vec![0, 1]);
    assert_eq!(call.genotype.format.pl_as_i32(), vec![45, 3, 0]);
    assert_eq!(call.genotype.format.gq.as_i32(), 3);
    assert!(
        !call.annotation_likelihoods.is_empty(),
        "per-variant annotation object must be attached"
    );

    let unique: BTreeSet<usize> = call
        .annotation_likelihoods
        .iter()
        .map(|c| c.read_index.get())
        .collect();
    kv("annotation_unique_n", unique.len().to_string());
    assert_eq!(unique.len(), 1, "not overlap-of-six; Java n=1");
    let idx = *unique.iter().next().expect("idx");
    let rec = outcome
        .genotyping_reads
        .get(idx)
        .expect("genotyping read")
        .as_ref();
    let qn = String::from_utf8_lossy(rec.qname()).into_owned();
    kv(
        "annotation_read",
        format!(
            "QNAME={qn} FLAG={} MAPQ={} strand={}",
            rec.flags(),
            rec.mapq(),
            if rec.flags() & FLAG_REVERSE != 0 {
                "rev"
            } else {
                "fwd"
            }
        ),
    );
    assert_eq!(qn, JAVA_READ);
    assert_eq!(rec.flags(), 83);
    assert_eq!(rec.mapq(), 44);

    let var_end = TARGET + 2;
    let overlap_n = outcome
        .genotyping_reads
        .iter()
        .filter(|r| java_alignment_read_overlaps_interval(r, TARGET, var_end, 2))
        .count();
    kv("overlap_n", overlap_n.to_string());
    assert_eq!(overlap_n, 6, "overlap-of-six still exists and is not used");
    assert_ne!(unique.len(), overlap_n);

    let region_dp = coverage_evidence_count(
        &outcome.genotyping_reads,
        &outcome.read_likelihoods,
        TARGET,
        var_end,
        2,
    );
    kv("region_wide_coverage_n", region_dp.to_string());
    let stored_unique: BTreeSet<usize> = outcome
        .read_likelihoods
        .iter()
        .map(|c| c.read_index.get())
        .collect();
    kv("region_wide_unique_n", stored_unique.len().to_string());
    assert_eq!(
        stored_unique.len(),
        2,
        "6R.185 stored poorly-modeled survivors"
    );
    assert_eq!(
        region_dp, 1,
        "6R.185 stored covering this deletion is 1 of 2 survivors"
    );

    let site_dp = coverage_evidence_count(
        &outcome.genotyping_reads,
        &call.annotation_likelihoods,
        TARGET,
        var_end,
        2,
    );
    let site_mqs =
        rms_mapping_quality_sample_mapqs(&outcome.genotyping_reads, &call.annotation_likelihoods);
    let site_mq = rms_mapping_quality_raw(&site_mqs);
    kv("site_coverage_n", site_dp.to_string());
    kv("site_mq", format!("{site_mq:?}"));
    assert_eq!(site_dp, 1);
    assert_eq!(site_mqs, vec![44]);
    assert_eq!(site_mq.map(|t| t.2), Some(44.0));

    let (full_ref, full_pad) = outcome.assembly.event_map_reference();
    let apply_pad = outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.as_ref().map(|g| g.start.get()))
        .unwrap_or(full_pad);
    let evidence = HcStrandBiasLikelihoods {
        reads: &outcome.genotyping_reads,
        likelihoods: &call.annotation_likelihoods,
        haplotypes: &outcome.assembly.haplotypes,
        contig: &covering.contig,
        ref_bytes: outcome.assembly.reference_bases(),
        pad_start_1based: apply_pad,
        full_ref_bytes: full_ref,
        full_pad_1based: full_pad,
        max_mnp_distance: outcome.assembly.max_mnp_distance(),
        emit_spanning_dels: true,
    };
    let sor_table = strand_bias_contingency_table(
        &evidence,
        TARGET,
        MERGED_REF,
        MERGED_ALT,
        STRAND_ODDS_RATIO_MIN_COUNT,
    );
    kv(
        "site_sor_table",
        format!(
            "[{},{};{},{}]",
            sor_table.0, sor_table.1, sor_table.2, sor_table.3
        ),
    );
    assert_eq!(sor_table, (0, 0, 0, 1), "ALT_REV only");
    assert!((java_calculate_sor(0, 0, 0, 1) - 1.6094379124341003).abs() < 1e-12);

    let emitted =
        try_emit_call_region_variants(covering, &outcome, "NA12878", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let rec = emitted
        .iter()
        .find(|r| r.position == TARGET && r.reference == MERGED_REF)
        .expect("record");
    let sample = rec.samples.first().expect("sample");
    assert_eq!(
        sample.gt.as_ref().map(|g| g.to_string()).as_deref(),
        Some("1/1")
    );
    assert_eq!(sample.ad.as_deref(), Some(&[0u32, 1][..]));
    assert_eq!(sample.dp.map(|d| d as i32), Some(1));
    assert_eq!(sample.gq.map(|g| g as i32), Some(3));
    assert_eq!(sample.pl.as_deref(), Some(&[45u32, 3, 0][..]));
    assert!((rec.quality.unwrap_or(0.0) - 35.48).abs() < 0.05);
    assert_eq!(info_i32(&rec.info, "DP"), Some(1));
    let mq = info_f64(&rec.info, "MQ").unwrap_or(-1.0);
    let sor = info_f64(&rec.info, "SOR").unwrap_or(-1.0);
    kv("emit_mq", format!("{mq}"));
    kv("emit_sor", format!("{sor}"));
    assert!((mq - 44.0).abs() < 0.005, "MQ={mq}");
    assert!(
        (sor - 1.6094379124341003).abs() < 1e-3 || (sor - 1.609).abs() < 0.002,
        "SOR={sor}"
    );
    assert!(!info_has(&rec.info, "InbreedingCoeff"));

    let closed_specs = parse_intervals_cli_string(&dict, "2:92305500-92305850").expect("closed");
    let closed_walk = traverse_assembly_region_walker(
        &dict,
        &closed_specs,
        &ref_fasta,
        &bam,
        &ReadFilterParams::gatk_standard_hc(),
        &WalkerTraversalConfig::gatk_haplotype_caller_production(100),
    )
    .expect("closed walk");
    let closed_regions = flatten_assembly_regions(&closed_walk);
    let closed_covering = closed_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_SNP
                && r.end.get() >= CLOSED_SNP
        })
        .expect("closed covering");
    let closed_outcome = HaplotypeCallerEngine::call_region(
        closed_covering,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("closed call")
    .expect("closed outcome");
    let closed_emitted = try_emit_call_region_variants(
        closed_covering,
        &closed_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("closed emit");
    let closed_rec = closed_emitted
        .iter()
        .find(|r| r.position == CLOSED_SNP)
        .expect("closed G/T");
    assert_eq!(info_i32(&closed_rec.info, "DP"), Some(3), "6R.174 stays");
    assert!(
        !info_has(&closed_rec.info, "InbreedingCoeff"),
        "6R.178 stays"
    );
    let closed_sor = info_f64(&closed_rec.info, "SOR").unwrap_or(-1.0);
    assert!((closed_sor - 0.6931471805599453).abs() < 1e-6 || (closed_sor - 0.693).abs() < 0.002);
    kv(
        "closed_snp",
        format!(
            "DP={} SOR={closed_sor} IC_omitted",
            info_i32(&closed_rec.info, "DP").unwrap_or(-1)
        ),
    );
}
