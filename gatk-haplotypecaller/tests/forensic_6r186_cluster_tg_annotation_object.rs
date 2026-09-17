//! 6R.186: construct the Java-equivalent per-variant annotation `AlleleLikelihoods`
//! on the strict-Java cluster-TG early-template path (production).
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! Target `2:92307333 T/G`. Consumes the 6R.185 stored n=2 hap matrix.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r186_cluster_tg_annotation_object -- --nocapture --test-threads=1
//! HOLDOUT_6R186=1 cargo test -p gatk-haplotypecaller --test holdout_6r186_cluster_tg_annotation -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{build_per_haplotype_variation_events, VariationEvent};
use gatk_haplotypecaller::hc_allele_mapping::create_allele_mapper_with_events;
use gatk_haplotypecaller::hc_genotyping_engine::java_alignment_read_overlaps_interval;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, marginalize_rows_to_biallelic_alleles,
    region_likelihoods_to_rows, traverse_assembly_region_walker, try_emit_call_region_variants,
    AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition, HaplotypeCallerEngine,
    ReadFilterParams, RegionReadLikelihood, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use rust_htslib::bam::Record;
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92307200-92307550";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_307_333;
const MERGED_REF: &str = "T";
const MERGED_ALT: &str = "G";
const CLOSED_INDEL: u64 = 92_307_324;
const CLOSED_SNP: u64 = 92_305_634;
const CLOSED_MID: u64 = 92_316_347;
const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const JAVA_MQ44_QNAME: &str = "H06JUADXX130110:1:1101:10052:88682";
const JAVA_SURVIVOR2_QNAME: &str = "H06HDADXX130110:1:1101:10061:17286";
const MARGIN: i32 = 2;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R186\t{key}\t{}", value.as_ref());
}

fn qname(rec: &Record) -> String {
    String::from_utf8_lossy(rec.qname()).into_owned()
}

fn unique_likelihood_indices(likelihoods: &[RegionReadLikelihood]) -> BTreeSet<usize> {
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
fn forensic_6r186_source_contracts() {
    kv("java_pin", JAVA_PIN);
    kv(
        "production_change",
        "cluster-TG early-template constructs annotation_likelihoods via marginalize then retainEvidence",
    );
    let early = include_str!("../src/hc_genotyping_engine/genotype_site_early_template.rs");
    let fin = include_str!("../src/hc_genotyping_engine/genotype_finalize.rs");
    let pipe = include_str!("../src/hc_genotyping_engine/genotype_site_pipeline.rs");
    let ge = include_str!("../src/hc_genotyping_engine/mod.rs");
    let emit = include_str!("../src/region_vcf_emit.rs");
    let ann = include_str!("../src/variant_site_hc_annotations.rs");
    let engine = include_str!("../src/engine.rs");
    let ll = include_str!("../src/engine_likelihoods.rs");

    assert!(
        early.contains("fn annotation_likelihoods_from_stored_haplotypes")
            && early.contains("marginalize_rows_to_biallelic_alleles")
            && early.contains("likelihood_subset_for_event"),
        "helper must run Java-order marginalize then retainEvidence"
    );
    let helper = early
        .split("fn annotation_likelihoods_from_stored_haplotypes")
        .nth(1)
        .expect("helper");
    assert!(
        !helper.contains("92307333")
            && !helper.contains(JAVA_MQ44_QNAME)
            && !helper.contains(JAVA_SURVIVOR2_QNAME),
        "helper has no coordinate/QNAME special case"
    );
    let tg_arm = early
        .split("if is_cluster_tg_snp(&event)")
        .nth(1)
        .expect("tg arm");
    let tg_body = tg_arm
        .split("if is_cluster_tc_snp(&event)")
        .next()
        .expect("tg body");
    assert!(
        tg_body.contains("annotation_likelihoods_from_stored_haplotypes")
            && tg_body.contains("with_annotation_likelihoods"),
        "cluster-TG early-template attaches the constructed object"
    );
    let fin_fn = fin
        .split("fn finish_strict_java_shaped_site_call")
        .nth(1)
        .expect("finish fn");
    let fin_body = fin_fn
        .split("fn finalize_strict_java_variation_genotype")
        .next()
        .expect("body");
    assert!(
        !fin_body.contains("with_annotation_likelihoods"),
        "shared finish helper must not globally attach annotation_likelihoods"
    );
    assert!(
        pipe.contains("with_annotation_likelihoods")
            && ge.contains("fn per_variant_annotation_likelihoods"),
        "6R.180 SiteScore path stays"
    );
    assert!(
        emit.contains("if call.annotation_likelihoods.is_empty()")
            && emit.contains("likelihoods: &call.annotation_likelihoods"),
        "emit still consumes the attached object when present"
    );
    assert!(
        engine.contains("fn apply_java_order_normalize_and_filter")
            && engine.contains("refresh_region_read_likelihoods"),
        "6R.185 stored-matrix lifecycle stays closed"
    );
    assert!(
        ann.contains("6R.174: INFO DP is Java `Coverage.evidenceCount`")
            && ann.contains("mq_rms_of_sample_evidence")
            && ann.contains("strand_bias_sample_counts"),
        "DP/MQ/SOR formulas must stay closed"
    );
    assert!(
        ll.contains("score_pairhmm_from_records_java_mate_contig"),
        "6R.176 mate-contig gate must stay closed"
    );
    assert!(
        !ann.contains("92307333") && !emit.contains("92307333") && !engine.contains("92307333"),
        "no locus-specific 6R.186 production patch"
    );
}

#[test]
fn forensic_6r186_cluster_tg_annotation_object() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
        return;
    }

    kv("java_pin", JAVA_PIN);
    kv("target", "2:92307333 T/G");
    kv(
        "production_change",
        "construct per-variant annotation AlleleLikelihoods on cluster-TG early-template",
    );
    kv("classification", "B — MISSING OBJECT CONSTRUCTION");

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

    let stored = unique_likelihood_indices(&outcome.read_likelihoods);
    kv("stored_hap_n", stored.len().to_string());
    assert_eq!(stored.len(), 2, "6R.185 stored hap matrix stays Java n=2");

    let (full_ref, full_pad) = outcome.assembly.event_map_reference();
    let apply_pad = outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.as_ref().map(|g| g.start.get()))
        .unwrap_or(full_pad);
    let hap_cache = build_per_haplotype_variation_events(
        &outcome.assembly.haplotypes,
        full_ref,
        full_pad,
        outcome.assembly.max_mnp_distance(),
        &covering.contig,
    );
    let mapping = create_allele_mapper_with_events(
        &VariationEvent::from_alleles(&covering.contig, TARGET, MERGED_REF, MERGED_ALT),
        TARGET,
        &outcome.assembly.haplotypes,
        apply_pad,
        outcome.assembly.reference_bases(),
        outcome.assembly.max_mnp_distance(),
        true,
        Some(&hap_cache),
    );
    let rows =
        region_likelihoods_to_rows(&outcome.read_likelihoods, outcome.assembly.haplotypes.len());
    let marg = marginalize_rows_to_biallelic_alleles(
        &rows,
        &mapping.ref_haplotype_indices,
        &mapping.alt_haplotype_indices,
    );
    let marg_unique: BTreeSet<usize> = marg.iter().map(|r| r.read_index).collect();
    kv("after_marginalize_n", marg_unique.len().to_string());
    assert_eq!(
        marg_unique.len(),
        2,
        "marginalize constructs a new n=2 allele object"
    );
    assert!(
        marg.iter()
            .all(|r| r.haplotype_log10_likelihoods.len() == 2),
        "marginalized columns are REF/ALT"
    );

    let retain: BTreeSet<usize> = marg_unique
        .iter()
        .copied()
        .filter(|&idx| {
            outcome
                .genotyping_reads
                .get(idx)
                .is_some_and(|r| java_alignment_read_overlaps_interval(r, TARGET, TARGET, MARGIN))
        })
        .collect();
    kv("after_retain_evidence_n", retain.len().to_string());
    assert_eq!(
        retain.len(),
        1,
        "Java-equivalent retainEvidence on the marginalized object is n=1"
    );

    let call = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(TARGET)
                && c.event.ref_allele == MERGED_REF
                && c.event.alt_allele == MERGED_ALT
        })
        .expect("T/G");
    assert_eq!(call.genotype.format.ad_as_i32(), vec![0, 1]);
    assert_eq!(call.genotype.format.dp.as_i32(), 1);
    assert_eq!(call.genotype.format.gq.as_i32(), 3);
    assert_eq!(call.genotype.format.pl_as_i32(), vec![45, 3, 0]);

    let attached = unique_likelihood_indices(&call.annotation_likelihoods);
    kv("annotation_likelihoods_n", attached.len().to_string());
    assert!(
        !call.annotation_likelihoods.is_empty(),
        "annotation_likelihoods must be present"
    );
    assert_eq!(attached.len(), 1, "attached object unique n=1");
    assert_eq!(
        attached, retain,
        "attached membership is retainEvidence n=1"
    );

    let keep_idx = *attached.iter().next().expect("n=1");
    let rec = outcome.genotyping_reads.get(keep_idx).expect("read");
    kv(
        "retained_read",
        format!(
            "{} FLAG={} MAPQ={} start={}",
            qname(rec),
            rec.flags(),
            rec.mapq(),
            rec.pos() + 1
        ),
    );
    assert_eq!(qname(rec), JAVA_MQ44_QNAME);
    assert_eq!(rec.flags(), 83);
    assert_eq!(rec.mapq(), 44);
    assert_eq!(rec.pos() + 1, 92_307_292);

    let marg_row = marg
        .iter()
        .find(|r| r.read_index == keep_idx)
        .expect("marg row");
    let ll_ref = marg_row
        .haplotype_log10_likelihoods
        .first()
        .copied()
        .unwrap_or(f64::NEG_INFINITY);
    let ll_alt = marg_row
        .haplotype_log10_likelihoods
        .get(1)
        .copied()
        .unwrap_or(f64::NEG_INFINITY);
    let best = if ll_alt > ll_ref { "G" } else { "T" };
    kv(
        "retained_best_allele",
        format!("{best} T*={ll_ref:.6} G={ll_alt:.6}"),
    );
    assert_eq!(best, "G");

    let dropped: Vec<_> = stored
        .difference(&attached)
        .filter_map(|&idx| outcome.genotyping_reads.get(idx))
        .collect();
    assert_eq!(dropped.len(), 1);
    assert_eq!(qname(dropped[0]), JAVA_SURVIVOR2_QNAME);
    assert_eq!(dropped[0].flags(), 83);
    assert_eq!(dropped[0].mapq(), 40);
    kv(
        "dropped_by_retain",
        format!(
            "{} FLAG={} MAPQ={} start={}",
            qname(dropped[0]),
            dropped[0].flags(),
            dropped[0].mapq(),
            dropped[0].pos() + 1
        ),
    );

    let indel = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(CLOSED_INDEL)
                && c.event.ref_allele == "TTC"
                && c.event.alt_allele == "T"
        })
        .expect("TTC/T");
    assert_eq!(
        unique_likelihood_indices(&indel.annotation_likelihoods).len(),
        1,
        "6R.180 TTC/T per-variant path stays"
    );

    let emitted =
        try_emit_call_region_variants(covering, &outcome, "NA12878", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let rec = emitted
        .iter()
        .find(|r| {
            r.position == TARGET
                && r.reference == MERGED_REF
                && r.alternate.first().map(String::as_str) == Some(MERGED_ALT)
        })
        .expect("T/G record");
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
    let fs = info_f64(&rec.info, "FS").unwrap_or(-1.0);
    kv("target_info", format!("DP=1 MQ={mq} SOR={sor} FS={fs}"));
    assert!((mq - 44.0).abs() < 0.005, "MQ={mq}");
    assert!(
        (sor - 1.6094379124341003).abs() < 1e-3 || (sor - 1.609).abs() < 0.002,
        "SOR={sor}"
    );

    let closed_indel = emitted
        .iter()
        .find(|r| r.position == CLOSED_INDEL && r.reference == "TTC")
        .expect("TTC/T");
    assert_eq!(info_i32(&closed_indel.info, "DP"), Some(1));
    let closed_mq = info_f64(&closed_indel.info, "MQ").unwrap_or(-1.0);
    let closed_sor = info_f64(&closed_indel.info, "SOR").unwrap_or(-1.0);
    assert!((closed_mq - 44.0).abs() < 0.005);
    assert!((closed_sor - 1.6094379124341003).abs() < 1e-3 || (closed_sor - 1.609).abs() < 0.002);

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
    let closed_sample = closed_rec.samples.first().expect("sample");
    assert_eq!(closed_sample.dp.map(|d| d as i32), Some(2));
    assert_eq!(info_i32(&closed_rec.info, "DP"), Some(3));
    let closed_snp_sor = info_f64(&closed_rec.info, "SOR").unwrap_or(-1.0);
    assert!((closed_snp_sor - 0.693).abs() < 0.002);
    assert!(!info_has(&closed_rec.info, "InbreedingCoeff"));

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
    let mid_covering = mid_regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= CLOSED_MID
                && r.end.get() >= CLOSED_MID
        })
        .expect("mid covering");
    let mid_outcome = HaplotypeCallerEngine::call_region(
        mid_covering,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("mid call")
    .expect("mid outcome");
    let mid_emitted = try_emit_call_region_variants(
        mid_covering,
        &mid_outcome,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("mid emit");
    let mid_rec = mid_emitted
        .iter()
        .find(|r| r.position == CLOSED_MID && r.reference == "G")
        .expect("2:92316347 G/A");
    let mid_fs = info_f64(&mid_rec.info, "FS").unwrap_or(-1.0);
    let mid_sor = info_f64(&mid_rec.info, "SOR").unwrap_or(-1.0);
    let mid_mq = info_f64(&mid_rec.info, "MQ").unwrap_or(-1.0);
    kv(
        "closed_92316347",
        format!("FS={mid_fs} SOR={mid_sor} MQ={mid_mq}"),
    );
    assert!((mid_fs - 0.0).abs() < 1e-9, "FS={mid_fs}");
    assert!((mid_sor - 1.179).abs() < 0.002, "SOR={mid_sor}");
    assert!((mid_mq - 40.25).abs() < 0.005, "MQ={mid_mq}");
}
