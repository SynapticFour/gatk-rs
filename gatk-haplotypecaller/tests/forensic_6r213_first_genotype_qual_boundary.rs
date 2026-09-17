//! 6R.213: first causal genotype/QUAL/QD arrow at `20:29455015 G/T`.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE. The first proven split is the genotype-likelihood
//! object for a true biallelic `G/T` site. QUAL/GQ/QD are downstream of PL.
//! Neighbor `20:29455019 G/A` already matches and is not this merge.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r213_first_genotype_qual_boundary -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{
    build_per_haplotype_variation_events, merged_alleles_for_genotyping,
    variation_events_at_position_from_cache,
};
use gatk_haplotypecaller::hc_allele_mapping::create_allele_mapper_with_events;
use gatk_haplotypecaller::variant_site_hc_annotations::qual_from_af_calculation;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, HcGenotypingConfig, ReadFilterParams, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_455_015;
const SIB: u64 = 29_455_019;
const TARGET_REF: &str = "G";
const TARGET_ALT: &str = "T";

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R213\t{key}\t{}", value.as_ref());
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

fn fmt_gls(gls: &[f64]) -> String {
    gls.iter()
        .map(|x| format!("{x:.8}"))
        .collect::<Vec<_>>()
        .join(",")
}

fn fmt_gt(s: &gatk_core::io::vcf::SampleData) -> String {
    s.gt.as_ref()
        .map(|g| {
            let sep = if g.phased { "|" } else { "/" };
            g.alleles
                .iter()
                .map(|a| a.to_string())
                .collect::<Vec<_>>()
                .join(sep)
        })
        .unwrap_or_default()
}

#[test]
fn forensic_6r213_java_record_pin() {
    kv("java_pin", JAVA_PIN);
    kv(
        "java_vcf",
        "20:29455015 G/T QUAL=61.64 FILTER=. AC=1 AF=0.500 AN=2 BaseQRankSum=1.820 DP=65 ExcessHet=0.0000 FS=0.000 MLEAC=1 MLEAF=0.500 MQ=43.37 MQRankSum=-3.327 QD=0.96 ReadPosRankSum=0.635 SOR=0.748 GT:AD:DP:GQ:PL 0/1:57,7:64:69:69,0,2140",
    );
    kv(
        "java_sib_29455019",
        "20:29455019 G/A QUAL=637.64 QD=10.28 GT:AD:DP:GQ:PL 0/1:39,23:62:99:645,0,1147",
    );
    kv("java_pl", "69,0,2140");
    kv("java_gq", "69");
    kv("java_gt", "0/1");
    kv("java_ad", "57,7");
    kv("java_fmt_dp", "64");
    kv("java_qual", "61.64");
    kv("java_qd", "0.96");
    kv("java_info_dp", "65");
}

#[test]
fn forensic_6r213_qual_from_java_pl_reproduces_printed_qual() {
    let java_gls = [-6.9_f64, 0.0, -214.0];
    let q = qual_from_af_calculation(&java_gls).expect("qual");
    assert!(
        (q - 61.64).abs() < 0.01,
        "AF calc on Java PL GLs must reproduce QUAL 61.64, got {q}"
    );
    assert!((q / 64.0 - 0.96).abs() < 0.01);
}

#[test]
fn forensic_6r213_qual_from_rust_pl_reproduces_current_qual() {
    let rust_gls = [-12.2_f64, 0.0, -230.4];
    let q = qual_from_af_calculation(&rust_gls).expect("qual");
    assert!(
        (q - 114.64).abs() < 0.05,
        "AF calc on Rust PL GLs must reproduce QUAL 114.64, got {q}"
    );
}

#[test]
fn forensic_6r213_live_first_causal_arrow() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }
    kv("java_pin", JAVA_PIN);
    kv(
        "production_change",
        "NONE in 6R.213; closed upstream by 6R.218 GRAPH_STATE_DIVERGENCE",
    );
    kv(
        "classification",
        "CLOSED_AFTER_6R218 (PL/GQ/QUAL/QD now match Java; first split was cyclic RT extract)",
    );
    kv("first_divergent_vcf_field", "none at 20:29455015");

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
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .collect();
    assert_eq!(covering.len(), 1, "one ActiveFull covers the target");
    let region = covering[0];

    let outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");

    let (full_ref, full_pad) = outcome.assembly.event_map_reference();
    let hap_events = build_per_haplotype_variation_events(
        &outcome.assembly.haplotypes,
        full_ref,
        full_pad,
        outcome.assembly.max_mnp_distance(),
        &region.contig,
    );
    let emit_spanning = !HcGenotypingConfig::strict_java().disable_spanning_event_genotyping;
    let cache_only = variation_events_at_position_from_cache(&hap_events, TARGET, emit_spanning);
    kv(
        "cache_only_events",
        cache_only
            .iter()
            .map(|e| format!("{}/{}", e.ref_allele, e.alt_allele))
            .collect::<Vec<_>>()
            .join(";"),
    );
    assert_eq!(cache_only.len(), 1);
    assert_eq!(cache_only[0].ref_allele, "G");
    assert_eq!(cache_only[0].alt_allele, "T");
    assert_eq!(
        merged_alleles_for_genotyping(&cache_only, TARGET),
        None,
        "production EventMap at 29455015 is biallelic G/T, not merged G/A,T"
    );
    let sib = variation_events_at_position_from_cache(&hap_events, SIB, emit_spanning);
    kv(
        "sib_29455019_events",
        sib.iter()
            .map(|e| format!("{}/{}", e.ref_allele, e.alt_allele))
            .collect::<Vec<_>>()
            .join(";"),
    );

    let apply_bases = outcome.assembly.apply_bases_shared();
    let apply_pad = outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.map(|g| g.start_1based()))
        .unwrap_or(full_pad);
    let mapping = create_allele_mapper_with_events(
        &cache_only[0],
        TARGET,
        &outcome.assembly.haplotypes,
        apply_pad,
        apply_bases.as_ref(),
        outcome.assembly.max_mnp_distance(),
        emit_spanning,
        Some(&hap_events),
    );
    kv(
        "mapper_ref_n",
        mapping.ref_haplotype_indices.len().to_string(),
    );
    kv(
        "mapper_alt_n",
        mapping.alt_haplotype_indices.len().to_string(),
    );
    kv("hap_n", outcome.assembly.haplotypes.len().to_string());
    kv(
        "read_likelihood_n",
        outcome.read_likelihoods.len().to_string(),
    );
    kv(
        "genotyping_read_n",
        outcome.genotyping_reads.len().to_string(),
    );

    let site = outcome
        .genotyped_calls
        .iter()
        .find(|c| {
            c.event.start_1based.get() == TARGET
                && c.event.ref_allele == TARGET_REF
                && c.event.alt_allele == TARGET_ALT
        })
        .expect("G/T site call");
    kv(
        "call_alleles",
        format!("{}/{}", site.event.ref_allele, site.event.alt_allele),
    );
    kv(
        "extra_alt_alleles",
        if site.extra_alt_alleles.is_empty() {
            "(none)".to_string()
        } else {
            site.extra_alt_alleles.join(",")
        },
    );
    kv(
        "post_merge_unused_alt_subset",
        site.post_merge_unused_alt_subset.to_string(),
    );
    assert!(!site.post_merge_unused_alt_subset);
    assert!(site.extra_alt_alleles.is_empty());
    kv(
        "raw_gls",
        fmt_gls(&site.genotype.genotype_log10_likelihoods),
    );
    let pl = site.genotype.format.pl_as_i32();
    kv(
        "pl",
        pl.iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(","),
    );
    kv("gq", site.genotype.format.gq.as_i32().to_string());
    kv(
        "ad",
        site.genotype
            .format
            .ad_as_i32()
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(","),
    );
    kv("fmt_dp", site.genotype.format.dp.as_i32().to_string());
    kv(
        "qual_log10_p_error",
        match site.qual_log10_p_error {
            Some(x) => format!("{x:.12}"),
            None => "None".to_string(),
        },
    );
    assert_eq!(pl, vec![69, 0, 2140], "6R.218 restored Java PL");
    assert_eq!(site.genotype.format.gq.as_i32(), 69);
    assert_eq!(site.genotype.format.ad_as_i32(), vec![57, 7]);
    assert_eq!(site.genotype.format.dp.as_i32(), 64);

    let gls = &site.genotype.genotype_log10_likelihoods;
    let q_from_gls = qual_from_af_calculation(gls).expect("qual from emitted GLs");
    kv("qual_from_emitted_gls", format!("{q_from_gls:.8}"));
    kv(
        "qd_from_emitted_gls_over_64",
        format!("{:.8}", q_from_gls / 64.0),
    );
    assert!(
        (q_from_gls - 62.11).abs() < 0.05,
        "AF calc on raw emitted GLs is ~62.11; VCF QUAL 61.64 is the PL round-trip, got {q_from_gls}"
    );
    let java_gls = [-6.9_f64, 0.0, -214.0];
    let q_java = qual_from_af_calculation(&java_gls).expect("java qual");
    kv("qual_from_java_pl_gls", format!("{q_java:.8}"));
    assert!((q_java - 61.64).abs() < 0.01);

    let emitted =
        try_emit_call_region_variants(region, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .unwrap_or_default();
    let sib_vcf = emitted
        .iter()
        .find(|r| r.position == SIB && r.reference == "G" && r.alternate.iter().any(|a| a == "A"))
        .expect("emitted G/A sibling");
    assert_eq!(
        sib_vcf.samples[0].pl.as_deref(),
        Some([645, 0, 1147].as_slice())
    );
    assert!((sib_vcf.quality.expect("sib QUAL") - 637.64).abs() < 0.01);
    let vcf = emitted
        .iter()
        .find(|r| {
            r.position == TARGET
                && r.reference == TARGET_REF
                && r.alternate.iter().any(|a| a == TARGET_ALT)
        })
        .expect("emitted G/T");
    kv("vcf_qual", format!("{}", vcf.quality.unwrap_or(-1.0)));
    let s = &vcf.samples[0];
    kv("vcf_gt", fmt_gt(s));
    kv(
        "vcf_ad",
        s.ad.as_ref()
            .map(|xs| {
                xs.iter()
                    .map(|x| x.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            })
            .unwrap_or_default(),
    );
    kv("vcf_dp", s.dp.map(|x| x.to_string()).unwrap_or_default());
    kv("vcf_gq", s.gq.map(|x| format!("{x}")).unwrap_or_default());
    kv(
        "vcf_pl",
        s.pl.as_ref()
            .map(|xs| {
                xs.iter()
                    .map(|x| x.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            })
            .unwrap_or_default(),
    );
    kv(
        "info_dp",
        info_i32(&vcf.info, "DP")
            .map(|x| x.to_string())
            .unwrap_or_default(),
    );
    kv(
        "info_qd",
        info_f64(&vcf.info, "QD")
            .map(|x| format!("{x}"))
            .unwrap_or_default(),
    );
    kv(
        "info_mq",
        info_f64(&vcf.info, "MQ")
            .map(|x| format!("{x:.2}"))
            .unwrap_or_default(),
    );
    kv(
        "info_fs",
        info_f64(&vcf.info, "FS")
            .map(|x| format!("{x}"))
            .unwrap_or_default(),
    );
    kv(
        "info_sor",
        info_f64(&vcf.info, "SOR")
            .map(|x| format!("{x}"))
            .unwrap_or_default(),
    );
    kv(
        "info_baseq",
        info_f64(&vcf.info, "BaseQRankSum")
            .map(|x| format!("{x}"))
            .unwrap_or_default(),
    );
    kv(
        "info_mqrs",
        info_f64(&vcf.info, "MQRankSum")
            .map(|x| format!("{x}"))
            .unwrap_or_default(),
    );
    kv(
        "info_rprs",
        info_f64(&vcf.info, "ReadPosRankSum")
            .map(|x| format!("{x}"))
            .unwrap_or_default(),
    );

    assert_eq!(fmt_gt(s), "0/1");
    assert_eq!(s.pl.as_deref(), Some([69, 0, 2140].as_slice()));
    assert_eq!(s.gq.map(|g| g as i32), Some(69));
    assert_eq!(s.ad.as_deref(), Some([57, 7].as_slice()));
    assert_eq!(s.dp, Some(64));
    assert!((vcf.quality.expect("QUAL") - 61.64).abs() < 0.05);
    assert_eq!(info_i32(&vcf.info, "DP"), Some(65));
    let qd = info_f64(&vcf.info, "QD").expect("QD");
    assert!((qd - 0.96).abs() < 0.02, "QD follows QUAL/64, got {qd}");
    assert!((info_f64(&vcf.info, "MQ").expect("MQ") - 43.37).abs() < 0.005);
    assert!(info_f64(&vcf.info, "FS").expect("FS").abs() < 0.01);
    assert!((info_f64(&vcf.info, "SOR").expect("SOR") - 0.748).abs() < 0.01);
    assert!((info_f64(&vcf.info, "BaseQRankSum").expect("BQ") - 1.820).abs() < 0.01);
    assert!((info_f64(&vcf.info, "MQRankSum").expect("MQRS") + 3.327).abs() < 0.01);
    assert!((info_f64(&vcf.info, "ReadPosRankSum").expect("RPRS") - 0.635).abs() < 0.01);
    kv(
        "causal_proof",
        "6R.218 cycle abort: PL/GQ/QUAL/QD now match Java 69,0,2140 / 69 / 61.64 / 0.96; GT/AD/DP/INFO unchanged",
    );
}
