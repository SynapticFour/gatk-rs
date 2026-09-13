//! 6R.185: restore Java-order normalize + poorly-modeled filter after the
//! final P12-cluster raw refresh (production).
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r185_post_refresh_java_order_filter -- --nocapture --test-threads=1
//! HOLDOUT_6R185=1 cargo test -p gatk-haplotypecaller --test holdout_6r185_post_refresh_filter -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_haplotypecaller::{
    begin_likelihood_pipeline_observe, begin_poorly_modeled_observe, call_disposition,
    flatten_assembly_regions, take_likelihood_pipeline_cells, take_likelihood_pipeline_snaps,
    take_poorly_modeled_observe, traverse_assembly_region_walker, try_emit_call_region_variants,
    AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition, HaplotypeCallerEngine,
    LikelihoodPipelineCell, ReadFilterParams, RegionReadLikelihood, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
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
const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const JAVA_MQ44_QNAME: &str = "H06JUADXX130110:1:1101:10052:88682";
const JAVA_SURVIVOR2_QNAME: &str = "H06HDADXX130110:1:1101:10061:17286";

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R185\t{key}\t{}", value.as_ref());
}

fn qname(rec: &Record) -> String {
    String::from_utf8_lossy(rec.qname()).into_owned()
}

fn unique_likelihood_indices(likelihoods: &[RegionReadLikelihood]) -> BTreeSet<usize> {
    likelihoods.iter().map(|c| c.read_index.get()).collect()
}

fn unique_from_cells(cells: &[LikelihoodPipelineCell], seq: u32) -> BTreeSet<usize> {
    cells
        .iter()
        .filter(|c| c.seq == seq)
        .map(|c| c.read_index)
        .collect()
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
fn forensic_6r185_source_reuses_existing_filter_after_last_p12_refresh() {
    kv("java_pin", JAVA_PIN);
    let engine = include_str!("../src/engine.rs");
    let early = include_str!("../src/hc_genotyping_engine/genotype_site_early_template.rs");
    let fin = include_str!("../src/hc_genotyping_engine/genotype_finalize.rs");
    let emit = include_str!("../src/region_vcf_emit.rs");
    let ann = include_str!("../src/variant_site_hc_annotations.rs");

    assert!(
        engine.contains("fn apply_java_order_normalize_and_filter")
            && engine.contains("normalize_region_read_likelihoods(read_likelihoods, &norm_haps)")
            && engine.contains("filter_normalized_region_read_likelihoods"),
        "production reuses existing normalize + poorly-modeled filter"
    );
    assert!(
        engine.contains("fn refresh_region_read_likelihoods")
            && engine.matches("refresh_region_read_likelihoods(").count() >= 2,
        "P12 refresh itself remains"
    );

    let outcome_move = engine
        .find("Ok(Some(CallRegionOutcome {")
        .expect("outcome move");
    let last_refresh = engine[..outcome_move]
        .rfind("read_likelihoods = ll;")
        .expect("last refresh assignment");
    let after_last = &engine[last_refresh..outcome_move];
    assert!(
        after_last.contains("apply_java_order_normalize_and_filter"),
        "Java-order normalize+filter runs after the last P12 refresh assignment"
    );
    assert!(
        after_last.contains("strict_java_p12_cluster_span"),
        "post-refresh repair is gated on the P12 cluster path, not a locus"
    );
    assert!(
        !after_last.contains("annotation_likelihoods"),
        "this arrow does not construct annotation_likelihoods"
    );

    let helper = engine
        .split("fn apply_java_order_normalize_and_filter")
        .nth(1)
        .expect("helper");
    assert!(
        !helper.contains("92307333") && !helper.contains(JAVA_MQ44_QNAME),
        "helper has no coordinate/QNAME special case"
    );
    assert!(
        !ann.contains("92307333") && !emit.contains("92307333") && !engine.contains("92307333"),
        "no locus-specific production patch"
    );
    assert!(
        early.contains("if is_cluster_tg_snp(&event)")
            && fin.contains("fn finish_strict_java_shaped_site_call"),
        "cluster-TG path still does not attach annotation here"
    );
}

#[test]
fn forensic_6r185_stored_membership_after_final_p12_filter() {
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
        "restore Java-order normalize+filter after last P12 refresh",
    );

    let dict = gatk_core::reference::SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs =
        gatk_core::reference::parse_intervals_cli_string(&dict, INTERVAL).expect("interval");
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

    begin_poorly_modeled_observe();
    begin_likelihood_pipeline_observe();
    let outcome = HaplotypeCallerEngine::call_region(
        covering,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("outcome");
    let observed = take_poorly_modeled_observe();
    let snaps = take_likelihood_pipeline_snaps();
    let cells = take_likelihood_pipeline_cells();

    for s in &snaps {
        kv(
            "LIFE_SNAP",
            format!(
                "seq={} stage={} n_reads={} n_haps={} unique={}",
                s.seq,
                s.stage,
                s.n_reads,
                s.n_haps,
                unique_from_cells(&cells, s.seq).len()
            ),
        );
    }

    let last = snaps.last().expect("last snap");
    kv("last_snap_stage", last.stage);
    assert_eq!(
        last.stage, "filter",
        "authoritative matrix is the post-refresh poorly-modeled filter"
    );
    let last_unique = unique_from_cells(&cells, last.seq);
    kv("last_filter_unique_n", last_unique.len().to_string());

    let stored = unique_likelihood_indices(&outcome.read_likelihoods);
    kv("stored_unique_n", stored.len().to_string());
    kv(
        "stored_hap_n",
        outcome.assembly.haplotypes.len().to_string(),
    );
    assert_eq!(stored.len(), 2, "stored genotyping evidence is Java n=2");
    assert_eq!(last_unique.len(), 2);
    assert_eq!(last_unique, stored);

    let last_pass = observed.iter().map(|r| r.pass).max().unwrap_or(0);
    let last_rows: Vec<_> = observed.iter().filter(|r| r.pass == last_pass).collect();
    let rust_keep_n = last_rows.iter().filter(|r| r.rust_keep).count();
    let extra_n = last_rows.iter().filter(|r| r.extra_retain).count();
    kv(
        "last_filter_pass",
        format!("pass={last_pass} rust_keep={rust_keep_n} extra={extra_n}"),
    );
    assert_eq!(rust_keep_n, 2);
    assert_eq!(extra_n, 0);

    let mut keep_ids = Vec::new();
    for &idx in &stored {
        let rec = outcome.genotyping_reads.get(idx).expect("read");
        keep_ids.push((qname(rec), rec.flags(), rec.mapq(), rec.pos() + 1));
        kv(
            "STORED_KEEP",
            format!(
                "{} FLAG={} MAPQ={} start={}",
                qname(rec),
                rec.flags(),
                rec.mapq(),
                rec.pos() + 1
            ),
        );
    }
    assert!(keep_ids
        .iter()
        .any(|(q, f, m, s)| { q == JAVA_MQ44_QNAME && *f == 83 && *m == 44 && *s == 92_307_292 }));
    assert!(keep_ids.iter().any(|(q, f, m, s)| {
        q == JAVA_SURVIVOR2_QNAME && *f == 83 && *m == 40 && *s == 92_307_338
    }));

    let call = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based == GenomePosition::new_1based(TARGET)
                && c.event.ref_allele == MERGED_REF
                && c.event.alt_allele == MERGED_ALT
        })
        .expect("T/G");
    assert!(
        call.annotation_likelihoods.is_empty(),
        "6R.181 B is not this arrow: annotation_likelihoods stays empty"
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

    let dp = info_i32(&rec.info, "DP");
    let mq = info_f64(&rec.info, "MQ");
    let sor = info_f64(&rec.info, "SOR");
    let fs = info_f64(&rec.info, "FS");
    kv(
        "target_downstream",
        format!(
            "stored_n=2 annotation_likelihoods=empty FORMAT=GT:AD:DP:GQ:PL 1/1:0,1:1:3:45,3,0 QUAL={:.2} INFO DP={:?} MQ={:?} SOR={:?} FS={:?}",
            rec.quality.unwrap_or(0.0),
            dp,
            mq,
            sor,
            fs
        ),
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

    let closed_specs =
        gatk_core::reference::parse_intervals_cli_string(&dict, "2:92305500-92305850")
            .expect("closed");
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
}
