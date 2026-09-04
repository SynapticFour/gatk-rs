//! 6R.74 coordinate-free: PairHMM input contract after quality preparation.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//! `PairHMMLikelihoodCalculationEngine.modifyReadQualities`
//! → `StandardPairHMMInputScoreImputator.impute`
//! → `PairHMM.computeLog10Likelihoods`.
//!
//! First remaining input field after 6R.73: BAM `BI`/`BD` must survive
//! `finalizeRegion` clipping. Java `ClippingOp.applyHardClipBases` copies and
//! slices those tags; Rust `replace_record_body` previously dropped them, so
//! production PairHMM saw Q45.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r74_pairhmm_input_contract
//! ```

use gatk_haplotypecaller::indel_gop_from_optional_tag;
use gatk_haplotypecaller::pairhmm_log10::GATK_PARITY_DEFAULT_GCP;
use gatk_haplotypecaller::pairhmm_qual::MIN_USABLE_Q_SCORE;
use gatk_haplotypecaller::pcr_error_model::{apply_pcr_error_model, PcrErrorModel};
use gatk_haplotypecaller::read_unclip::{
    hard_clip_low_qual_ends, revert_soft_clipped_bases, HcSoftclipPolicy,
};
use gatk_haplotypecaller::{prepare_read_quals_for_pairhmm_inplace, HcLikelihoodEngineConfig};
use rust_htslib::bam;
use rust_htslib::bam::record::{Aux, Cigar, CigarString};

const ANTI_SEQ: &[u8] = b"ACGT";
const ANTI_BI: u8 = 30;
const ANTI_BD: u8 = 25;
const ANTI_BQ: u8 = 31;
const LOW_TAIL_BQ: u8 = 2;
const MAPQ: u8 = 60;
const MIN_TAIL: u8 = 9;

#[derive(Clone, Debug, Eq, PartialEq)]
struct PairHmmInput {
    read_bases: Vec<u8>,
    hap_bases: Vec<u8>,
    bq: Vec<u8>,
    iq: Vec<u8>,
    dq: Vec<u8>,
    gcp: Vec<u8>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum InputField {
    ReadBases,
    HapBases,
    Bq,
    Iq,
    Dq,
    Gcp,
    ReadLen,
    HapLen,
}

fn first_differing_field(java: &PairHmmInput, rust: &PairHmmInput) -> Option<InputField> {
    if java.read_bases != rust.read_bases {
        return Some(if java.read_bases.len() != rust.read_bases.len() {
            InputField::ReadLen
        } else {
            InputField::ReadBases
        });
    }
    if java.hap_bases != rust.hap_bases {
        return Some(if java.hap_bases.len() != rust.hap_bases.len() {
            InputField::HapLen
        } else {
            InputField::HapBases
        });
    }
    if java.bq != rust.bq {
        return Some(InputField::Bq);
    }
    if java.iq != rust.iq {
        return Some(InputField::Iq);
    }
    if java.dq != rust.dq {
        return Some(InputField::Dq);
    }
    if java.gcp != rust.gcp {
        return Some(InputField::Gcp);
    }
    None
}

fn phred_to_fastq(q: &[u8]) -> String {
    q.iter()
        .map(|&b| char::from(b.saturating_add(33)))
        .collect()
}

fn aux_phred(rec: &bam::Record, tag: &[u8]) -> Option<Vec<u8>> {
    match rec.aux(tag) {
        Ok(Aux::String(s)) => Some(s.bytes().map(|b| b.saturating_sub(33)).collect()),
        _ => None,
    }
}

fn make_read(seq: &[u8], bq: &[u8], cigar: CigarString, bi: &[u8], bd: &[u8]) -> bam::Record {
    let mut rec = bam::Record::new();
    rec.set(b"r1", Some(&cigar), seq, bq);
    rec.set_mapq(MAPQ);
    rec.set_pos(100);
    rec.set_tid(0);
    rec.push_aux(b"BI", Aux::String(&phred_to_fastq(bi)))
        .expect("BI");
    rec.push_aux(b"BD", Aux::String(&phred_to_fastq(bd)))
        .expect("BD");
    rec
}

/// Java `modifyReadQualities` + `StandardPairHMMInputScoreImputator` on already-clipped bases.
fn java_pairhmm_input(
    read_bases: &[u8],
    raw_bq: &[u8],
    mapq: u8,
    bi: Option<&[u8]>,
    bd: Option<&[u8]>,
    hap_bases: &[u8],
) -> PairHmmInput {
    let cfg = HcLikelihoodEngineConfig::default();
    let mut bq = raw_bq.to_vec();
    prepare_read_quals_for_pairhmm_inplace(&mut bq, mapq, &cfg);
    let mut iq = indel_gop_from_optional_tag(bi, read_bases.len()).unwrap();
    let mut dq = indel_gop_from_optional_tag(bd, read_bases.len()).unwrap();
    apply_pcr_error_model(read_bases, &mut iq, &mut dq, PcrErrorModel::Conservative);
    for q in iq.iter_mut().chain(dq.iter_mut()) {
        if *q < MIN_USABLE_Q_SCORE {
            *q = MIN_USABLE_Q_SCORE;
        }
    }
    PairHmmInput {
        read_bases: read_bases.to_vec(),
        hap_bases: hap_bases.to_vec(),
        bq,
        iq,
        dq,
        gcp: vec![GATK_PARITY_DEFAULT_GCP; read_bases.len()],
    }
}

/// Production `score_read_against_haplotypes` prep, using tags on the *scored* record.
fn rust_pairhmm_input_from_record(rec: &bam::Record, hap_bases: &[u8]) -> PairHmmInput {
    let cfg = HcLikelihoodEngineConfig::default();
    let read_bases = rec.seq().as_bytes();
    let mut bq = rec.qual().to_vec();
    prepare_read_quals_for_pairhmm_inplace(&mut bq, rec.mapq(), &cfg);
    let bi = aux_phred(rec, b"BI");
    let bd = aux_phred(rec, b"BD");
    let mut iq = indel_gop_from_optional_tag(bi.as_deref(), read_bases.len()).unwrap();
    let mut dq = indel_gop_from_optional_tag(bd.as_deref(), read_bases.len()).unwrap();
    apply_pcr_error_model(&read_bases, &mut iq, &mut dq, PcrErrorModel::Conservative);
    PairHmmInput {
        read_bases,
        hap_bases: hap_bases.to_vec(),
        bq,
        iq,
        dq,
        gcp: vec![GATK_PARITY_DEFAULT_GCP; rec.seq_len()],
    }
}

/// Legacy 6R.73-on-finalized-record: missing BI/BD → Q45.
fn rust_legacy_q45_after_dropped_tags(rec: &bam::Record, hap_bases: &[u8]) -> PairHmmInput {
    let cfg = HcLikelihoodEngineConfig::default();
    let read_bases = rec.seq().as_bytes();
    let mut bq = rec.qual().to_vec();
    prepare_read_quals_for_pairhmm_inplace(&mut bq, rec.mapq(), &cfg);
    let mut iq = indel_gop_from_optional_tag(None, read_bases.len()).unwrap();
    let mut dq = indel_gop_from_optional_tag(None, read_bases.len()).unwrap();
    apply_pcr_error_model(&read_bases, &mut iq, &mut dq, PcrErrorModel::Conservative);
    PairHmmInput {
        read_bases,
        hap_bases: hap_bases.to_vec(),
        bq,
        iq,
        dq,
        gcp: vec![GATK_PARITY_DEFAULT_GCP; rec.seq_len()],
    }
}

fn hap_ref() -> Vec<u8> {
    ANTI_SEQ[..3].to_vec()
}

fn hap_alt() -> Vec<u8> {
    let mut h = hap_ref();
    h[1] = b'T';
    h
}

#[test]
fn first_differing_field_is_iq_when_bi_is_dropped() {
    let rec = make_read(
        ANTI_SEQ,
        &[ANTI_BQ, ANTI_BQ, ANTI_BQ, LOW_TAIL_BQ],
        CigarString::from(vec![Cigar::Match(4)]),
        &[ANTI_BI; 4],
        &[ANTI_BD; 4],
    );
    let clipped = hard_clip_low_qual_ends(&rec, MIN_TAIL);
    assert_eq!(clipped.seq().as_bytes(), b"ACG");
    let hap = hap_ref();
    let java = java_pairhmm_input(
        &clipped.seq().as_bytes(),
        clipped.qual(),
        MAPQ,
        Some(&[ANTI_BI, ANTI_BI, ANTI_BI]),
        Some(&[ANTI_BD, ANTI_BD, ANTI_BD]),
        &hap,
    );
    let legacy = rust_legacy_q45_after_dropped_tags(&clipped, &hap);
    assert_eq!(
        first_differing_field(&java, &legacy),
        Some(InputField::Iq),
        "anti-masking: dropped BI → Q45 IQ vs sliced BI"
    );
    assert_ne!(java.iq[0], legacy.iq[0]);
    assert_eq!(java.iq[0], ANTI_BI);
    assert_eq!(legacy.iq[0], 40, "Q45 PCR-capped to cache[1]");
}

#[test]
fn hard_clip_pairhmm_inputs_match_java_after_bi_bd_preserve() {
    let rec = make_read(
        ANTI_SEQ,
        &[ANTI_BQ, ANTI_BQ, ANTI_BQ, LOW_TAIL_BQ],
        CigarString::from(vec![Cigar::Match(4)]),
        &[ANTI_BI; 4],
        &[ANTI_BD; 4],
    );
    let clipped = hard_clip_low_qual_ends(&rec, MIN_TAIL);
    let hap = hap_alt();
    let java = java_pairhmm_input(
        &clipped.seq().as_bytes(),
        clipped.qual(),
        MAPQ,
        aux_phred(&clipped, b"BI").as_deref(),
        aux_phred(&clipped, b"BD").as_deref(),
        &hap,
    );
    let rust = rust_pairhmm_input_from_record(&clipped, &hap);
    assert_eq!(first_differing_field(&java, &rust), None);
    assert_eq!(rust.iq, vec![ANTI_BI, ANTI_BI, ANTI_BI]);
    assert_eq!(rust.dq, vec![ANTI_BD, ANTI_BD, ANTI_BD]);
    assert_eq!(rust.gcp, vec![10, 10, 10]);
    assert_eq!(rust.read_bases, b"ACG");
    assert_eq!(rust.hap_bases, hap);
}

#[test]
fn revert_softclip_keeps_full_bi_bd_as_pairhmm_iq_dq() {
    let rec = make_read(
        ANTI_SEQ,
        &[ANTI_BQ; 4],
        CigarString::from(vec![Cigar::SoftClip(1), Cigar::Match(3)]),
        &[ANTI_BI; 4],
        &[ANTI_BD; 4],
    );
    let reverted = revert_soft_clipped_bases(&rec);
    let hap = hap_ref();
    let rust = rust_pairhmm_input_from_record(&reverted, &hap);
    let java = java_pairhmm_input(
        ANTI_SEQ,
        &[ANTI_BQ; 4],
        MAPQ,
        Some(&[ANTI_BI; 4]),
        Some(&[ANTI_BD; 4]),
        &hap,
    );
    assert_eq!(first_differing_field(&java, &rust), None);
    assert_eq!(rust.iq[0], ANTI_BI);
    assert_eq!(rust.dq[0], ANTI_BD);
}

#[test]
fn pair_ordering_is_read_outer_haplotype_inner() {
    let haps = [hap_ref(), hap_alt()];
    let mut pairs = Vec::new();
    for (ri, _) in [0, 1].iter().enumerate() {
        for (hi, h) in haps.iter().enumerate() {
            pairs.push((ri, hi, h.as_slice()));
        }
    }
    assert_eq!(
        pairs.iter().map(|(r, h, _)| (*r, *h)).collect::<Vec<_>>(),
        vec![(0, 0), (0, 1), (1, 0), (1, 1)]
    );
}

#[test]
fn hc_default_keeps_softclips_before_pairhmm() {
    let policy = HcSoftclipPolicy::haplotype_caller_defaults();
    assert!(!policy.dont_use_soft_clipped_bases);
}

#[test]
fn iq_dq_floor_is_noop_when_gop_at_least_six() {
    assert!(ANTI_BI > MIN_USABLE_Q_SCORE);
    assert!(ANTI_BD > MIN_USABLE_Q_SCORE);
}

/// Live capture for the report. Not a coordinate-only proof.
///
/// ```text
/// HOLDOUT_6R74=1 cargo test -p gatk-haplotypecaller --test forensic_6r74_pairhmm_input_contract canonical_live -- --nocapture
/// ```
#[test]
fn canonical_live_pairhmm_matrix_dimensions() {
    if std::env::var("HOLDOUT_6R74").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R74=1");
        return;
    }
    use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
    use gatk_haplotypecaller::assembly_region_finalize::{
        finalize_region_reads_for_assembly, gatk_min_tail_quality_for_assembly,
    };
    use gatk_haplotypecaller::{
        call_disposition, flatten_assembly_regions, take_colocated_merge_numerics,
        traverse_assembly_region_walker, try_emit_call_region_variants,
        AssemblyRegionCallDisposition, CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams,
        WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
    };
    use std::collections::HashSet;
    use std::path::Path;

    const INTERVAL: &str = "20:29455000-29456500";
    const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
    const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
    const POS_SNP: u64 = 29_456_344;

    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: live BAM/ref missing");
        return;
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
        .expect("lifecycle T/C");
    let n_reads: HashSet<usize> = outcome
        .read_likelihoods
        .iter()
        .map(|c| c.read_index.get())
        .collect();
    let n_haps: HashSet<usize> = outcome
        .read_likelihoods
        .iter()
        .map(|c| c.haplotype_index.get())
        .collect();
    let ref_len = outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .map(|h| h.bases.len())
        .unwrap_or(0);
    let n_long_alt = outcome
        .assembly
        .haplotypes
        .iter()
        .filter(|h| !h.is_reference && h.bases.len() > ref_len.saturating_add(8))
        .count();
    let n_bi_on_geno = outcome
        .genotyping_reads
        .iter()
        .filter(|r| aux_phred(r, b"BI").is_some())
        .count();
    let finalized = finalize_region_reads_for_assembly(
        &covering[0].reads,
        covering[0],
        true,
        gatk_min_tail_quality_for_assembly(10),
        false,
    );
    let n_final = finalized.len();
    let mut n_bi_final = 0usize;
    let mut n_bi_len_ok = 0usize;
    let mut n_iq_ne_q45 = 0usize;
    let mut first_iq_diff: Option<(usize, u8, u8)> = None;
    for rec in &finalized {
        let seq_len = rec.seq().as_bytes().len();
        if let Some(bi) = aux_phred(rec, b"BI") {
            n_bi_final += 1;
            if bi.len() == seq_len {
                n_bi_len_ok += 1;
                if let Some(i) = bi.iter().position(|&q| q != 45) {
                    n_iq_ne_q45 += 1;
                    if first_iq_diff.is_none() {
                        first_iq_diff = Some((i, 45, bi[i]));
                    }
                }
            }
        }
    }
    let live = take_colocated_merge_numerics();
    let numerics = live.iter().find(|n| n.loc == POS_SNP);
    eprintln!(
        "6R.74 live: PL={:?} AD={:?} QUAL={:?} geno_reads={} ll_reads={} ll_haps={} assembly_haps={} long_alt>ref+8={} cells={} overlap={:?} bi_on_geno={} finalized={} bi_final={} bi_len_ok={} iq_ne_q45={} first_iq_diff={:?}",
        vcf.samples[0].pl,
        vcf.samples[0].ad,
        vcf.quality,
        outcome.genotyping_reads.len(),
        n_reads.len(),
        n_haps.len(),
        outcome.assembly.haplotypes.len(),
        n_long_alt,
        outcome.read_likelihoods.len(),
        numerics.map(|n| (n.n_reads, n.pool_sizes.clone())),
        n_bi_on_geno,
        n_final,
        n_bi_final,
        n_bi_len_ok,
        n_iq_ne_q45,
        first_iq_diff,
    );
    assert_eq!(
        n_bi_final, n_final,
        "Java ClippingOp keeps BI through finalize"
    );
    assert_eq!(n_bi_len_ok, n_final, "clipped BI length equals seq length");
    assert_eq!(
        n_reads.len().saturating_mul(n_haps.len()),
        outcome.read_likelihoods.len()
    );
}
