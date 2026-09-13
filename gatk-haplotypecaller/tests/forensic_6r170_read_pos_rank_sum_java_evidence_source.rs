//! 6R.170: ReadPosRankSum evidence is `region.reads` pileup, not Java RankSumTest.
//!
//! Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`, GKL 0.8.8.
//! Target `2:92316347 G/A` is a live fixture. Production change: NONE.
//! FORMAT/QUAL/FS/SOR/MQ stay closed. InbreedingCoeff is out of scope.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r170_read_pos_rank_sum_java_evidence_source -- --test-threads=1 --nocapture
//! HOLDOUT_6R170=1 cargo test -p gatk-haplotypecaller --test holdout_6r170_read_pos_rank_sum -- --nocapture --test-threads=1
//! ```

use gatk_core::io::vcf::InfoValue;
use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::event_map::{build_per_haplotype_variation_events, VariationEvent};
use gatk_haplotypecaller::fragment_overlap::read_base_at_ref_coord_1based;
use gatk_haplotypecaller::hc_allele_mapping::create_allele_mapper_with_events;
use gatk_haplotypecaller::read_model::MAPPING_QUALITY_UNAVAILABLE;
use gatk_haplotypecaller::read_realignment::LOG_10_INFORMATIVE_THRESHOLD;
use gatk_haplotypecaller::variant_site_hc_annotations::HcStrandBiasLikelihoods;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, marginalize_rows_to_biallelic_alleles,
    region_likelihoods_to_rows, traverse_assembly_region_walker, try_emit_call_region_variants,
    AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition, HaplotypeCallerEngine,
    HcGenotypingConfig, ReadFilterParams, WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Record;
use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92316200-92316580";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_316_347;
const MERGED_REF: &str = "G";
const MERGED_ALT: &str = "A";
const FLAG_REVERSE: u16 = 0x10;
/// 6R.166 informative best-allele members (post-filter AlleleLikelihoods).
const JAVA_INFORMATIVE_MEMBERS: &[&str] = &[
    "H06HDADXX130110:1:1101:10034:45116 FLAG=99 strand=fwd allele=ALT",
    "H06HDADXX130110:2:1101:10025:49248 FLAG=99 strand=fwd allele=ALT",
    "H06HDADXX130110:2:1101:10046:78083 FLAG=147 strand=rev allele=ALT",
];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R170\t{key}\t{}", value.as_ref());
}

fn info_f64(info: &[InfoValue], key: &str) -> f64 {
    for v in info {
        if let InfoValue::Float(k, xs) = v {
            if k == key {
                return xs.first().copied().unwrap_or(0.0);
            }
        }
    }
    0.0
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

fn strand_label(rec: &Record) -> &'static str {
    if rec.flags() & FLAG_REVERSE != 0 {
        "rev"
    } else {
        "fwd"
    }
}

fn member_allele(qname: &str, rec: &Record, allele: &str) -> String {
    format!(
        "{qname} FLAG={} strand={} allele={allele}",
        rec.flags(),
        strand_label(rec)
    )
}

/// Java `RankSumTest.isUsableRead`: MQ != 0 and MQ != 255.
fn java_ranksum_mq_usable(rec: &Record) -> bool {
    let mq = rec.mapq();
    mq != 0 && mq != MAPPING_QUALITY_UNAVAILABLE
}

/// Java `ReadPosRankSumTest.getReadPosition`: `min(left, right)` including hard clips.
fn java_read_pos_rank_element(rec: &Record, vc_start_1based: i32) -> Option<f64> {
    let cig = rec.cigar();
    let mut leading_hard = 0i64;
    let mut trailing_hard = 0i64;
    if let Some(Cigar::HardClip(n)) = cig.iter().next() {
        leading_hard = i64::from(*n);
    }
    if let Some(Cigar::HardClip(n)) = cig.iter().last() {
        trailing_hard = i64::from(*n);
    }
    let vc = i64::from(vc_start_1based);
    let mut ref_pos = rec.pos() + 1;
    let mut q = 0usize;
    for c in cig.iter() {
        match c {
            Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                let n = i64::from(*n);
                if vc >= ref_pos && vc < ref_pos + n {
                    let idx = q + (vc - ref_pos) as usize;
                    let seq_len = rec.seq_len() as i64;
                    let left = leading_hard + idx as i64;
                    let right = seq_len - 1 - idx as i64 + trailing_hard;
                    return Some(left.min(right) as f64);
                }
                ref_pos += n;
                q += n as usize;
            }
            Cigar::Del(n) | Cigar::RefSkip(n) => {
                ref_pos += i64::from(*n);
            }
            Cigar::Ins(n) | Cigar::SoftClip(n) => {
                q += *n as usize;
            }
            Cigar::HardClip(_) | Cigar::Pad(_) => {}
        }
    }
    None
}

/// Coordinate-free: production still uses pileup offsets (6R.170 membership).
/// 6R.172 closed empty→0.0; this test does not pin that old sentinel.
#[test]
fn forensic_6r170_production_still_pileup_and_inserts_zero() {
    let ann = include_str!("../src/variant_site_hc_annotations.rs");
    assert!(
        ann.contains("read_offset_evidence_at_site")
            && ann.contains("read_pos_rank_sum::read_pos_rank_sum"),
        "ReadPosRankSum must still be computed from read_offset_evidence_at_site"
    );
    assert!(
        ann.contains("6R.166: FS/SOR use Java") && ann.contains("6R.167: MQ membership is Java"),
        "FS/SOR/MQ paths must stay closed"
    );
    let rp = include_str!("../src/annotator/plugins/read_pos_rank_sum.rs");
    let production_plugin = rp.split("#[cfg(test)]").next().expect("plugin body");
    assert!(
        production_plugin.contains("if alt.is_empty() || reference.is_empty()")
            && production_plugin.contains("return None"),
        "6R.172: empty lists are None; 6R.170 still documents pileup membership"
    );
    assert!(
        production_plugin.contains("result.z.is_nan()") && production_plugin.contains("None"),
        "6R.172: MannWhitneyU NaN is None"
    );
    let emit = include_str!("../src/region_vcf_emit.rs");
    assert!(
        emit.contains("if let Some(z) = ann.read_pos_rank_sum"),
        "6R.172: hc_info_values inserts ReadPosRankSum only for Some(z)"
    );
    assert!(
        !ann.contains("92316347") && !emit.contains("92316347") && !rp.contains("92316347"),
        "no locus-specific ReadPosRankSum patch"
    );
    let mw = include_str!("../src/mann_whitney_u.rs");
    assert!(
        mw.contains("if n1 == 0 || n2 == 0") && mw.contains("f64::NAN"),
        "MannWhitneyU.test still returns NaN on empty series (Java contract)"
    );
}

#[test]
fn forensic_6r170_read_pos_rank_sum_java_evidence_source() {
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing P12 ref/BAM");
        return;
    }

    kv(
        "java_pin",
        "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77",
    );
    kv("variant", format!("2:{TARGET} G/A"));
    kv(
        "production_change",
        "NONE (6R.170); emit predicate closed in 6R.172",
    );
    kv(
        "java_path",
        "StandardAnnotation → ReadPosRankSumTest → RankSumTest.fillQualsFromLikelihood → bestAllelesBreakingTies + isInformative + isUsableRead → MannWhitneyU NaN → emptyMap",
    );
    kv(
        "rust_path",
        "annotate_hc_variant_site → read_offset_evidence_at_site(region.reads pileup) → read_pos_rank_sum (6R.172 Option<f64>)",
    );

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
        .expect("ActiveFull covering target");

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
        .expect("genotyped G/A");
    let fmt = &call.genotype.format;
    assert_eq!(fmt.ad_as_i32(), vec![0, 3]);
    assert_eq!(fmt.pl_as_i32(), vec![135, 9, 0]);
    assert_eq!(fmt.dp.as_i32(), 3);
    assert_eq!(fmt.gq.as_i32(), 9);

    let (full_ref, full_pad) = outcome.assembly.event_map_reference();
    let apply_pad = outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.as_ref().map(|g| g.start.get()))
        .unwrap_or(full_pad);
    let cfg = HcGenotypingConfig::default();
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
    kv(
        "likelihood_object",
        "reuse outcome.read_likelihoods + genotyping_reads (post filterPoorlyModeledEvidence; same as 6R.166/167)",
    );

    let mut sample_ids = BTreeSet::new();
    let sample_qnames: BTreeSet<String> = {
        let mut seen = BTreeSet::new();
        let mut names = BTreeSet::new();
        for cell in evidence.likelihoods {
            let idx = cell.read_index.get();
            if !seen.insert(idx) {
                continue;
            }
            if let Some(rec) = evidence.reads.get(idx) {
                sample_ids.insert(format!(
                    "{} FLAG={}",
                    String::from_utf8_lossy(rec.qname()),
                    rec.flags()
                ));
                names.insert(format!(
                    "{} FLAG={} MAPQ={} CIGAR={} strand={}",
                    String::from_utf8_lossy(rec.qname()),
                    rec.flags(),
                    rec.mapq(),
                    rec.cigar(),
                    strand_label(rec)
                ));
            }
        }
        names
    };
    for row in &sample_qnames {
        kv("sample_evidence", row);
    }
    kv("sample_evidence_n", format!("{}", sample_qnames.len()));

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
    let rows = region_likelihoods_to_rows(evidence.likelihoods, evidence.haplotypes.len());
    let marg = marginalize_rows_to_biallelic_alleles(
        &rows,
        &mapping.ref_haplotype_indices,
        &mapping.alt_haplotype_indices,
    );

    let mut java_ref = Vec::new();
    let mut java_alt = Vec::new();
    let mut informative_members = Vec::new();
    for row in &marg {
        let Some(rec) = evidence.reads.get(row.read_index) else {
            continue;
        };
        let qname = String::from_utf8_lossy(rec.qname()).into_owned();
        let lls = &row.haplotype_log10_likelihoods;
        let ll_ref = lls.first().copied().unwrap_or(f64::NEG_INFINITY);
        let ll_alt = lls.get(1).copied().unwrap_or(f64::NEG_INFINITY);
        if !ll_ref.is_finite() && !ll_alt.is_finite() {
            kv(
                "likelihood_skip_nonfinite",
                format!("{qname} FLAG={} MAPQ={}", rec.flags(), rec.mapq()),
            );
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
        let usable = java_ranksum_mq_usable(rec);
        let pos = java_read_pos_rank_element(rec, TARGET as i32);
        kv(
            "likelihood_read",
            format!(
                "{qname} FLAG={} MAPQ={} strand={} CIGAR={} best={} informative={informative} mq_usable={usable} pos={pos:?} filterPoorlyModeledEvidence=survived gap={gap:.4}",
                rec.flags(),
                rec.mapq(),
                strand_label(rec),
                rec.cigar(),
                if best_is_ref { "REF" } else { "ALT" }
            ),
        );
        if !informative {
            continue;
        }
        informative_members.push(member_allele(
            &qname,
            rec,
            if best_is_ref { "REF" } else { "ALT" },
        ));
        if !usable {
            kv("ranksum_drop_mq", format!("{qname} MAPQ={}", rec.mapq()));
            continue;
        }
        let Some(p) = pos else {
            kv("ranksum_drop_pos", format!("{qname} no getReadPosition"));
            continue;
        };
        if best_is_ref {
            java_ref.push(p);
        } else {
            java_alt.push(p);
        }
    }
    assert_eq!(
        informative_members, JAVA_INFORMATIVE_MEMBERS,
        "Java RankSumTest fillQualsFromLikelihood informative membership must match 6R.166"
    );
    kv(
        "java_equivalent_counts",
        format!("REF={} ALT={}", java_ref.len(), java_alt.len()),
    );
    kv("java_equivalent_ref_pos", format!("{java_ref:?}"));
    kv("java_equivalent_alt_pos", format!("{java_alt:?}"));
    assert_eq!(
        java_ref.len(),
        0,
        "Java-equivalent REF RankSum evidence must be empty"
    );
    assert_eq!(
        java_alt.len(),
        3,
        "Java-equivalent ALT RankSum evidence must be 3"
    );

    let mut rust_ref = Vec::new();
    let mut rust_alt = Vec::new();
    let ref_b = MERGED_REF.as_bytes()[0];
    let alt_b = MERGED_ALT.as_bytes()[0];
    for rec in &covering.reads {
        let Some(base) = read_base_at_ref_coord_1based(rec, TARGET as i32) else {
            continue;
        };
        let offset = (TARGET as i64 - (rec.pos() + 1)) as f64;
        let allele = if base.eq_ignore_ascii_case(&alt_b) {
            rust_alt.push(offset);
            "ALT"
        } else if base.eq_ignore_ascii_case(&ref_b) {
            rust_ref.push(offset);
            "REF"
        } else {
            continue;
        };
        kv(
            "rust_pileup_read",
            format!(
                "{} FLAG={} MAPQ={} strand={} CIGAR={} allele={allele} pileup_offset={offset} in_sample_evidence={}",
                String::from_utf8_lossy(rec.qname()),
                rec.flags(),
                rec.mapq(),
                strand_label(rec),
                rec.cigar(),
                sample_ids.contains(&format!(
                    "{} FLAG={}",
                    String::from_utf8_lossy(rec.qname()),
                    rec.flags()
                ))
            ),
        );
    }
    kv(
        "rust_pileup_counts",
        format!("REF={} ALT={}", rust_ref.len(), rust_alt.len()),
    );
    assert_eq!(
        rust_ref.len(),
        2,
        "current Rust pileup still has 2 REF bases"
    );
    assert_eq!(
        rust_alt.len(),
        5,
        "current Rust pileup still has 5 ALT bases"
    );
    assert_ne!(
        (rust_ref.len(), rust_alt.len()),
        (java_ref.len(), java_alt.len()),
        "first arrow is evidence membership, not Mann-Whitney math"
    );

    // 6R.172 closed empty→0.0. Membership (pileup vs informative) remains the 6R.170 arrow.
    let plugin = include_str!("../src/annotator/plugins/read_pos_rank_sum.rs");
    let production_plugin = plugin.split("#[cfg(test)]").next().expect("plugin body");
    assert!(
        production_plugin.contains("if alt.is_empty() || reference.is_empty()")
            && production_plugin.contains("return None"),
        "Java-equivalent empty REF is None after 6R.172"
    );
    kv(
        "downstream_emit_predicate",
        "6R.172: Java-equivalent REF=0 ALT=3 → None → no insert. Live pileup REF=2 ALT=5 is a separate membership arrow.",
    );

    let emitted =
        try_emit_call_region_variants(covering, &outcome, "SAMPLE", DEFAULT_STAND_EMIT_CONFIDENCE)
            .expect("emit");
    let rec = emitted
        .iter()
        .find(|r| {
            r.position == TARGET
                && r.reference == MERGED_REF
                && r.alternate.first().map(String::as_str) == Some(MERGED_ALT)
        })
        .expect("emitted G/A");
    let fs = info_f64(&rec.info, "FS");
    let sor = info_f64(&rec.info, "SOR");
    let mq = info_f64(&rec.info, "MQ");
    let rp = info_f64(&rec.info, "ReadPosRankSum");
    kv(
        "rust_final",
        format!(
            "FS={fs} SOR={sor} MQ={mq} ReadPosRankSum={rp} has_key={} QUAL={:?}",
            info_has(&rec.info, "ReadPosRankSum"),
            rec.quality
        ),
    );
    assert!(info_has(&rec.info, "ReadPosRankSum"));
    assert!(
        (rp - 0.0).abs() < 1e-12,
        "production still emits 0, got {rp}"
    );
    assert!(fs < 0.02, "FS must stay Java 0");
    assert!((sor - 1.179).abs() < 0.001, "SOR must stay Java 1.179");
    assert!((mq - 40.25).abs() < 1e-12, "MQ must stay 40.25");
    assert!((rec.quality.unwrap_or(0.0) - 121.84).abs() < 0.02);
    assert!(
        !info_has(&rec.info, "InbreedingCoeff"),
        "6R.178: n=1 omits InbreedingCoeff"
    );
}
