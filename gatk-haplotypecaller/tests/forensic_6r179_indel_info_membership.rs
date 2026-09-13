//! 6R.179: first INFO causal arrow at `2:92307324 TTC/T` (proof-only).
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE.
//!
//! After 6R.178, this is the first FORMAT/QUAL-matched INFO numeric split.
//! Do not assume DP vs MQ vs SOR; dump live membership for all three.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r179_indel_info_membership -- --nocapture --test-threads=1
//! HOLDOUT_6R179=1 cargo test -p gatk-haplotypecaller --test holdout_6r179_indel_info -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{build_per_haplotype_variation_events, VariationEvent};
use gatk_haplotypecaller::hc_allele_mapping::create_allele_mapper_with_events;
use gatk_haplotypecaller::hc_genotyping_engine::{
    java_alignment_read_covers_variant_base, java_alignment_read_overlaps_interval,
};
use gatk_haplotypecaller::read_assembly_filter::{passes_assembly_read, AssemblyReadFilterConfig};
use gatk_haplotypecaller::read_model::MAPPING_QUALITY_UNAVAILABLE;
use gatk_haplotypecaller::read_realignment::LOG_10_INFORMATIVE_THRESHOLD;
use gatk_haplotypecaller::variant_site_hc_annotations::{
    coverage_evidence_count, rms_mapping_quality_raw, rms_mapping_quality_sample_mapqs,
    rms_mapping_quality_sample_reads, strand_bias_contingency_table, HcStrandBiasLikelihoods,
    STRAND_ODDS_RATIO_MIN_COUNT,
};
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, marginalize_rows_to_biallelic_alleles,
    region_likelihoods_to_rows, traverse_assembly_region_walker, try_emit_call_region_variants,
    AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition, HaplotypeCallerEngine,
    HcGenotypingConfig, ReadFilterParams, RegionReadLikelihood, WalkerTraversalConfig,
    DEFAULT_STAND_EMIT_CONFIDENCE,
};
use rust_htslib::bam::Record;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92307200-92307550";
const CLOSED_SNP: u64 = 92_305_634;
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_VCF_REL: &str = "parity/reports/6r43/p12_indel_mix/java.vcf";
const TARGET: u64 = 92_307_324;
const MERGED_REF: &str = "TTC";
const MERGED_ALT: &str = "T";
const FLAG_REVERSE: u16 = 0x10;
const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const MARGIN: i32 = 2;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R179\t{key}\t{}", value.as_ref());
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

fn qname(rec: &Record) -> String {
    String::from_utf8_lossy(rec.qname()).into_owned()
}

fn strand_label(rec: &Record) -> &'static str {
    if rec.flags() & FLAG_REVERSE != 0 {
        "rev"
    } else {
        "fwd"
    }
}

fn parse_info_map(info: &str) -> BTreeMap<String, String> {
    let mut out = BTreeMap::new();
    for piece in info.split(';') {
        if piece.is_empty() {
            continue;
        }
        match piece.split_once('=') {
            Some((k, v)) => {
                out.insert(k.to_string(), v.to_string());
            }
            None => {
                out.insert(piece.to_string(), String::new());
            }
        }
    }
    out
}

fn java_target_record(path: &Path) -> (BTreeMap<String, String>, String, String) {
    let text = std::fs::read_to_string(path).expect("java.vcf");
    for line in text.lines() {
        if line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 10 {
            continue;
        }
        if f[0] == "2" && f[1] == "92307324" && f[3] == MERGED_REF && f[4] == MERGED_ALT {
            return (parse_info_map(f[7]), f[8].to_string(), f[9].to_string());
        }
    }
    panic!("java.vcf missing 2:92307324 TTC/T");
}

fn mate_ok_java(rec: &Record) -> bool {
    if !rec.is_paired() || rec.is_mate_unmapped() || rec.is_unmapped() {
        return true;
    }
    rec.tid() == rec.mtid()
}

fn dump_read(prefix: &str, rec: &Record) {
    let cfg = AssemblyReadFilterConfig::gatk_defaults();
    eprintln!(
        "6R179\t{prefix}\tQNAME={} FLAG={} MAPQ={} strand={} start={} end={} CIGAR={} tid={} mtid={} mate_unmapped={} mate_ok={} asm_fail={}",
        qname(rec),
        rec.flags(),
        rec.mapq(),
        strand_label(rec),
        rec.pos() + 1,
        gatk_haplotypecaller::read_unclip::alignment_end_1based(rec),
        rec.cigar(),
        rec.tid(),
        rec.mtid(),
        rec.is_mate_unmapped(),
        mate_ok_java(rec),
        !passes_assembly_read(rec, &cfg),
    );
}

fn unique_likelihood_indices(likelihoods: &[RegionReadLikelihood]) -> BTreeSet<usize> {
    likelihoods.iter().map(|c| c.read_index.get()).collect()
}

fn real_pairhmm(likelihoods: &[RegionReadLikelihood], idx: usize) -> bool {
    likelihoods
        .iter()
        .any(|c| c.read_index.get() == idx && c.log10_likelihood != 0.0)
}

#[test]
fn forensic_6r179_source_contract_no_production_change() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    let ann = include_str!("../src/variant_site_hc_annotations.rs");
    let emit = include_str!("../src/region_vcf_emit.rs");
    let ll = include_str!("../src/engine_likelihoods.rs");
    assert!(
        ann.contains("6R.174: INFO DP is Java `Coverage.evidenceCount`")
            && ann.contains("coverage_evidence_count("),
        "6R.174 INFO DP path must stay closed"
    );
    assert!(
        ann.contains("6R.167: MQ membership is Java `sampleEvidence`")
            && ann.contains("mq_rms_of_sample_evidence"),
        "6R.167/168 MQ path must stay closed"
    );
    assert!(
        ann.contains("6R.166: FS/SOR use Java `getContingencyTable`")
            && ann.contains("strand_bias_sample_counts"),
        "6R.166 SOR path must stay closed"
    );
    assert!(
        ll.contains("score_pairhmm_from_records_java_mate_contig")
            && ll.contains("MATE_ON_SAME_CONTIG_OR_NO_MAPPED_MATE"),
        "6R.176 mate-contig PairHMM gate must stay closed"
    );
    assert!(
        emit.contains("INBREEDING_COEFF_MIN_SAMPLES")
            && emit.contains("n_genotypes >= INBREEDING_COEFF_MIN_SAMPLES"),
        "6R.178 InbreedingCoeff gate must stay closed"
    );
    assert!(
        !ann.contains("92307324") && !emit.contains("92307324") && !ll.contains("92307324"),
        "no locus-specific 6R.179 production patch"
    );
    let one_alt = java_calculate_sor(0, 0, 1, 0);
    kv("java_sor_table_0_0_1_0", format!("{one_alt:.6}"));
    assert!(
        (one_alt - 1.6094379124341003).abs() < 1e-12,
        "Java [0,0;1,0] prints SOR=1.609"
    );
}

#[test]
fn forensic_6r179_indel_info_membership() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    let java_vcf = root.join(JAVA_VCF_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
        return;
    }

    kv("java_pin", JAVA_PIN);
    kv("target", "2:92307324 TTC/T");
    kv("production_change", "NONE");

    if java_vcf.is_file() {
        let (jinfo, jfmt, jsample) = java_target_record(&java_vcf);
        kv("java_format_keys", &jfmt);
        kv("java_sample", jsample.clone());
        kv("java_info_dp", jinfo.get("DP").cloned().unwrap_or_default());
        kv("java_info_mq", jinfo.get("MQ").cloned().unwrap_or_default());
        kv(
            "java_info_sor",
            jinfo.get("SOR").cloned().unwrap_or_default(),
        );
        kv("java_info_fs", jinfo.get("FS").cloned().unwrap_or_default());
        kv(
            "java_rp",
            if jinfo.contains_key("ReadPosRankSum") {
                jinfo["ReadPosRankSum"].clone()
            } else {
                "OMITTED".into()
            },
        );
        kv(
            "java_ic",
            if jinfo.contains_key("InbreedingCoeff") {
                jinfo["InbreedingCoeff"].clone()
            } else {
                "OMITTED".into()
            },
        );
        assert_eq!(jinfo.get("DP").map(String::as_str), Some("1"));
        assert_eq!(jinfo.get("MQ").map(String::as_str), Some("44.00"));
        assert_eq!(jinfo.get("SOR").map(String::as_str), Some("1.609"));
        assert_eq!(jinfo.get("FS").map(String::as_str), Some("0.000"));
        assert!(!jinfo.contains_key("ReadPosRankSum"));
        assert!(!jinfo.contains_key("InbreedingCoeff"));
        let fields: Vec<&str> = jfmt.split(':').collect();
        let vals: Vec<&str> = jsample.split(':').collect();
        let get = |name: &str| {
            fields
                .iter()
                .zip(vals.iter())
                .find(|(k, _)| **k == name)
                .map(|(_, v)| *v)
                .unwrap_or("")
        };
        assert_eq!(get("GT"), "1/1");
        assert_eq!(get("AD"), "0,1");
        assert_eq!(get("DP"), "1");
        assert_eq!(get("GQ"), "3");
        assert_eq!(get("PL"), "45,3,0");
    }

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
    kv(
        "active_region",
        format!(
            "{}:{}-{}",
            covering.contig,
            covering.start.get(),
            covering.end.get()
        ),
    );

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
        .expect("genotyped TTC/T");
    let fmt = &call.genotype.format;
    assert_eq!(fmt.ad_as_i32(), vec![0, 1], "FORMAT AD closed");
    assert_eq!(fmt.pl_as_i32(), vec![45, 3, 0], "FORMAT PL closed");
    assert_eq!(fmt.dp.as_i32(), 1, "FORMAT DP closed");
    assert_eq!(fmt.gq.as_i32(), 3, "FORMAT GQ closed");

    let var_end = VariationEvent::vcf_end_1based(TARGET, MERGED_REF);
    kv("vcf_span", format!("{TARGET}-{var_end} (TTC deletion)"));
    kv("margin", MARGIN.to_string());

    let cfg = HcGenotypingConfig::default();
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
        likelihoods: &outcome.read_likelihoods,
        haplotypes: &outcome.assembly.haplotypes,
        contig: &covering.contig,
        ref_bytes: outcome.assembly.reference_bases(),
        pad_start_1based: apply_pad,
        full_ref_bytes: full_ref,
        full_pad_1based: full_pad,
        max_mnp_distance: outcome.assembly.max_mnp_distance(),
        emit_spanning_dels: !cfg.disable_spanning_event_genotyping,
    };

    let dp_n = coverage_evidence_count(
        evidence.reads,
        evidence.likelihoods,
        TARGET,
        var_end,
        cfg.informative_read_overlap_margin,
    );
    let mq_reads = rms_mapping_quality_sample_reads(evidence.reads, evidence.likelihoods);
    let mq_mqs = rms_mapping_quality_sample_mapqs(evidence.reads, evidence.likelihoods);
    let mq_raw = rms_mapping_quality_raw(&mq_mqs);
    let sor_table = strand_bias_contingency_table(
        &evidence,
        TARGET,
        MERGED_REF,
        MERGED_ALT,
        STRAND_ODDS_RATIO_MIN_COUNT,
    );
    kv("rust_coverage_n", dp_n.to_string());
    kv("rust_mq_sampleEvidence_n", mq_mqs.len().to_string());
    kv("rust_mq_mapqs", format!("{mq_mqs:?}"));
    kv("rust_mq_raw", format!("{mq_raw:?}"));
    kv(
        "rust_sor_table",
        format!(
            "[{},{};{},{}]",
            sor_table.0, sor_table.1, sor_table.2, sor_table.3
        ),
    );
    kv(
        "rust_sor_from_table",
        format!(
            "{:.6}",
            java_calculate_sor(sor_table.0, sor_table.1, sor_table.2, sor_table.3)
        ),
    );

    let ll_idx = unique_likelihood_indices(evidence.likelihoods);
    kv("region_reads", covering.reads.len().to_string());
    kv("genotyping_reads", evidence.reads.len().to_string());
    kv("likelihood_unique", ll_idx.len().to_string());

    let hap_cache = build_per_haplotype_variation_events(
        evidence.haplotypes,
        evidence.full_ref_bytes,
        evidence.full_pad_1based,
        evidence.max_mnp_distance,
        evidence.contig,
    );
    let mapping = create_allele_mapper_with_events(
        &VariationEvent::from_alleles(evidence.contig, TARGET, MERGED_REF, MERGED_ALT),
        TARGET,
        evidence.haplotypes,
        evidence.pad_start_1based,
        evidence.ref_bytes,
        evidence.max_mnp_distance,
        evidence.emit_spanning_dels,
        Some(&hap_cache),
    );
    kv(
        "haplotypes",
        format!(
            "n={} ref_pool={:?} alt_pool={:?}",
            evidence.haplotypes.len(),
            mapping
                .ref_haplotype_indices
                .iter()
                .map(|i| i.get())
                .collect::<Vec<_>>(),
            mapping
                .alt_haplotype_indices
                .iter()
                .map(|i| i.get())
                .collect::<Vec<_>>()
        ),
    );
    for (i, h) in evidence.haplotypes.iter().enumerate() {
        let events: Vec<String> = hap_cache
            .events_for(i)
            .iter()
            .map(|e| {
                format!(
                    "{}:{} {}/{}",
                    e.contig,
                    e.start_1based.get(),
                    e.ref_allele,
                    e.alt_allele
                )
            })
            .collect();
        kv(
            "hap",
            format!(
                "i={i} is_ref={} events={}",
                h.is_reference,
                events.join(",")
            ),
        );
    }

    let rows = region_likelihoods_to_rows(evidence.likelihoods, evidence.haplotypes.len());
    let marg = marginalize_rows_to_biallelic_alleles(
        &rows,
        &mapping.ref_haplotype_indices,
        &mapping.alt_haplotype_indices,
    );

    let orig_by_key: BTreeMap<(Vec<u8>, u16), &Record> = covering
        .reads
        .iter()
        .map(|r| ((r.qname().to_vec(), r.flags()), r.as_ref()))
        .collect();

    let mut overlap_m2 = BTreeSet::new();
    let mut overlap_m0 = BTreeSet::new();
    let mut covers_m2 = BTreeSet::new();
    let mut qname_n: BTreeMap<String, usize> = BTreeMap::new();
    for rec in covering.reads.iter() {
        dump_read("REGION", rec);
        *qname_n.entry(qname(rec)).or_insert(0) += 1;
    }
    for (name, n) in &qname_n {
        if *n > 1 {
            kv("multi_qname", format!("{name} n={n}"));
        }
    }

    for rec in evidence.reads {
        dump_read("GENOTYPING", rec);
    }

    for &idx in &ll_idx {
        let Some(rec) = evidence.reads.get(idx) else {
            continue;
        };
        let orig = orig_by_key
            .get(&(rec.qname().to_vec(), rec.flags()))
            .copied();
        let mate_ok = orig.map(mate_ok_java).unwrap_or_else(|| mate_ok_java(rec));
        let ov_m2 = java_alignment_read_overlaps_interval(rec, TARGET, var_end, MARGIN);
        let ov_m0 = java_alignment_read_overlaps_interval(rec, TARGET, var_end, 0);
        let cov_m2 = java_alignment_read_covers_variant_base(rec, TARGET, var_end, MARGIN);
        if ov_m2 {
            overlap_m2.insert(idx);
        }
        if ov_m0 {
            overlap_m0.insert(idx);
        }
        if cov_m2 {
            covers_m2.insert(idx);
        }
        kv(
            "LL_READ",
            format!(
                "idx={idx} QNAME={} FLAG={} MAPQ={} strand={} start={} end={} CIGAR={} mate_ok={} pairhmm_real={} overlap_m2={} overlap_m0={} covers_m2={} in_mq={}",
                qname(rec),
                rec.flags(),
                rec.mapq(),
                strand_label(rec),
                rec.pos() + 1,
                gatk_haplotypecaller::read_unclip::alignment_end_1based(rec),
                rec.cigar(),
                mate_ok,
                real_pairhmm(evidence.likelihoods, idx),
                ov_m2,
                ov_m0,
                cov_m2,
                mq_reads.iter().any(|r| r.qname() == rec.qname() && r.flags() == rec.flags()),
            ),
        );
    }

    kv("overlap_m2_n", overlap_m2.len().to_string());
    kv("overlap_m0_n", overlap_m0.len().to_string());
    kv("covers_m2_n", covers_m2.len().to_string());

    let overlap_mqs: Vec<u8> = overlap_m2
        .iter()
        .filter_map(|i| evidence.reads.get(*i))
        .map(|r| r.mapq())
        .filter(|&mq| mq != MAPPING_QUALITY_UNAVAILABLE)
        .collect();
    let covers_mqs: Vec<u8> = covers_m2
        .iter()
        .filter_map(|i| evidence.reads.get(*i))
        .map(|r| r.mapq())
        .filter(|&mq| mq != MAPPING_QUALITY_UNAVAILABLE)
        .collect();
    kv("overlap_m2_mapqs", format!("{overlap_mqs:?}"));
    kv("covers_m2_mapqs", format!("{covers_mqs:?}"));
    kv(
        "overlap_m2_rms",
        format!("{:?}", rms_mapping_quality_raw(&overlap_mqs)),
    );
    kv(
        "covers_m2_rms",
        format!("{:?}", rms_mapping_quality_raw(&covers_mqs)),
    );

    let mut full_inf = (0u32, 0u32, 0u32, 0u32);
    let mut overlap_inf = (0u32, 0u32, 0u32, 0u32);
    let mut covers_inf = (0u32, 0u32, 0u32, 0u32);
    for row in &marg {
        let Some(rec) = evidence.reads.get(row.read_index) else {
            continue;
        };
        let lls = &row.haplotype_log10_likelihoods;
        let ll_ref = lls.first().copied().unwrap_or(f64::NEG_INFINITY);
        let ll_alt = lls.get(1).copied().unwrap_or(f64::NEG_INFINITY);
        if !ll_ref.is_finite() && !ll_alt.is_finite() {
            continue;
        }
        let (best_is_ref, best, second) = if ll_ref > ll_alt {
            (true, ll_ref, ll_alt)
        } else if ll_alt > ll_ref {
            (false, ll_alt, ll_ref)
        } else {
            (true, ll_ref, ll_alt)
        };
        let gap = if second.is_finite() {
            best - second
        } else {
            f64::INFINITY
        };
        let informative = gap > LOG_10_INFORMATIVE_THRESHOLD;
        let allele = if !informative {
            "NON_INF"
        } else if best_is_ref {
            "REF"
        } else {
            "ALT"
        };
        kv(
            "MARG",
            format!(
                "idx={} QNAME={} FLAG={} MAPQ={} strand={} ll_ref={ll_ref:.3} ll_alt={ll_alt:.3} gap={gap:.3} allele={allele} overlap={} covers={}",
                row.read_index,
                qname(rec),
                rec.flags(),
                rec.mapq(),
                strand_label(rec),
                overlap_m2.contains(&row.read_index),
                covers_m2.contains(&row.read_index),
            ),
        );
        if !informative {
            continue;
        }
        let reverse = rec.flags() & FLAG_REVERSE != 0;
        let cell = match (best_is_ref, reverse) {
            (true, false) => (1, 0, 0, 0),
            (true, true) => (0, 1, 0, 0),
            (false, false) => (0, 0, 1, 0),
            (false, true) => (0, 0, 0, 1),
        };
        full_inf.0 += cell.0;
        full_inf.1 += cell.1;
        full_inf.2 += cell.2;
        full_inf.3 += cell.3;
        if overlap_m2.contains(&row.read_index) {
            overlap_inf.0 += cell.0;
            overlap_inf.1 += cell.1;
            overlap_inf.2 += cell.2;
            overlap_inf.3 += cell.3;
        }
        if covers_m2.contains(&row.read_index) {
            covers_inf.0 += cell.0;
            covers_inf.1 += cell.1;
            covers_inf.2 += cell.2;
            covers_inf.3 += cell.3;
        }
    }
    kv(
        "sor_full_region",
        format!(
            "{:?} -> {:.6}",
            full_inf,
            java_calculate_sor(full_inf.0, full_inf.1, full_inf.2, full_inf.3)
        ),
    );
    kv(
        "sor_overlap_m2",
        format!(
            "{:?} -> {:.6}",
            overlap_inf,
            java_calculate_sor(overlap_inf.0, overlap_inf.1, overlap_inf.2, overlap_inf.3)
        ),
    );
    kv(
        "sor_covers_m2",
        format!(
            "{:?} -> {:.6}",
            covers_inf,
            java_calculate_sor(covers_inf.0, covers_inf.1, covers_inf.2, covers_inf.3)
        ),
    );

    let emitted =
        try_emit_call_region_variants(covering, &outcome, "NA12878", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let rec = emitted
        .iter()
        .find(|r| r.position == TARGET && r.reference == MERGED_REF)
        .expect("emitted TTC/T");
    let sample = rec.samples.first().expect("sample");
    assert_eq!(
        sample.gt.as_ref().map(|g| g.to_string()).as_deref(),
        Some("1/1")
    );
    assert_eq!(sample.ad.as_deref(), Some(&[0u32, 1][..]));
    assert_eq!(sample.dp.map(|d| d as i32), Some(1));
    assert_eq!(sample.gq.map(|g| g as i32), Some(3));
    assert_eq!(sample.pl.as_deref(), Some(&[45u32, 3, 0][..]));
    assert_eq!(
        info_i32(&rec.info, "DP"),
        Some(1),
        "6R.180: emit INFO DP is per-variant n=1, not region coverage {dp_n}"
    );
    let mq = info_f64(&rec.info, "MQ").unwrap_or(-1.0);
    let sor = info_f64(&rec.info, "SOR").unwrap_or(-1.0);
    let fs = info_f64(&rec.info, "FS").unwrap_or(-1.0);
    kv(
        "rust_emit_dp",
        info_i32(&rec.info, "DP").unwrap_or(-1).to_string(),
    );
    kv("rust_emit_mq", format!("{mq}"));
    kv("rust_emit_sor", format!("{sor}"));
    kv("rust_emit_fs", format!("{fs}"));
    kv("rust_emit_qual", format!("{:?}", rec.quality));
    assert!(
        !info_has(&rec.info, "InbreedingCoeff"),
        "6R.178 stays closed"
    );
    assert!(fs < 0.02, "FS stays ~0");
    assert!(
        (rec.quality.unwrap_or(0.0) - 35.44).abs() < 0.05,
        "QUAL print-close"
    );

    // Closed SNP from 6R.174–6R.178 must not move.
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
    let closed_out = HaplotypeCallerEngine::call_region(
        closed_covering,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("closed call")
    .expect("closed outcome");
    let closed_emitted = try_emit_call_region_variants(
        closed_covering,
        &closed_out,
        "NA12878",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .expect("closed emit");
    let closed_rec = closed_emitted
        .iter()
        .find(|r| r.position == CLOSED_SNP)
        .expect("closed G/T");
    assert_eq!(info_i32(&closed_rec.info, "DP"), Some(3));
    assert!(!info_has(&closed_rec.info, "InbreedingCoeff"));
    let closed_sor = info_f64(&closed_rec.info, "SOR").unwrap_or(-1.0);
    assert!((closed_sor - 0.6931471805599453).abs() < 1e-6 || (closed_sor - 0.693).abs() < 0.002);

    assert_eq!(
        dp_n, 1,
        "6R.185 stored poorly-modeled survivors covering this deletion (n=1 of 2)"
    );
    assert_eq!(
        mq_mqs.len(),
        2,
        "6R.185 stored unique evidence n=2 (region-wide sampleEvidence)"
    );
    assert_eq!(
        sor_table,
        (0, 1, 0, 1),
        "region-wide stored after 6R.185; emit still uses 6R.180 per-variant object"
    );
    assert_eq!(overlap_m2.len(), 1);
    assert_eq!(
        overlap_m0.len(),
        1,
        "margin=0 still 1 of the 2 stored survivors: not a ±2 overlap-predicate miss"
    );
    assert_eq!(
        covers_m2.len(),
        1,
        "covering a query base still 1 of 2 stored survivors: not a deletion spanning-base miss"
    );
    let java_one = java_calculate_sor(0, 0, 0, 1);
    assert!((java_one - 1.6094379124341003).abs() < 1e-12);
    assert!(
        (java_calculate_sor(1, 2, 2, 2) - 0.36772478012531723).abs() < 1e-9,
        "Rust emitted SOR is calculateSOR on [1,2;2,2], not a formula split"
    );
    let one_mq = rms_mapping_quality_raw(&[44]);
    assert_eq!(one_mq.map(|t| t.2), Some(44.0));
    kv(
        "java_one_read",
        "H06JUADXX130110:1:1101:10052:88682 FLAG=83 MAPQ=44 ALT_REV reconstructs DP=1 MQ=44.00 SOR=1.609",
    );
    kv(
        "java_later_site_overlap6",
        "2:92307403 Java INFO DP=6 MQ=40.58 = RMS of the same 6 overlapping MAPQs; TTC/T n=1 is a per-variant subset, not 'only 6 reads exist'",
    );
    kv(
        "first_arrow",
        "6R.185 stored outcome.read_likelihoods is Java n=2 (overlap=1 at this deletion). TTC/T emit still uses the 6R.180 per-variant annotation object (n=1, MAPQ=44). Cluster-TG T/G still has empty annotation_likelihoods (later 6R.181 B).",
    );
}
