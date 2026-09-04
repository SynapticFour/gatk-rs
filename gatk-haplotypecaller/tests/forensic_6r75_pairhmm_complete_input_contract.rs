//! 6R.75 coordinate-free: complete PairHMM kernel-boundary input tuple.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//! `modifyReadQualities` → `StandardPairHMMInputScoreImputator.impute`
//! → `PairHMM.computeLog10Likelihoods` arrays, immediately before
//! `subComputeReadLikelihoodGivenHaplotypeLog10`.
//!
//! Rust: `score_read_against_haplotypes` scratch planes immediately before
//! Log10 / Logless / SIMD. Scoring is not entered as a proof.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r75_pairhmm_complete_input_contract
//! ```

use gatk_haplotypecaller::haplotype::{sort_haplotypes_assembly_result_order, Haplotype};
use gatk_haplotypecaller::indel_gop_from_optional_tag;
use gatk_haplotypecaller::pairhmm_log10::GATK_PARITY_DEFAULT_GCP;
use gatk_haplotypecaller::pairhmm_qual::MIN_USABLE_Q_SCORE;
use gatk_haplotypecaller::pcr_error_model::{apply_pcr_error_model, PcrErrorModel};
use gatk_haplotypecaller::read_unclip::hard_clip_low_qual_ends;
use gatk_haplotypecaller::{prepare_read_quals_for_pairhmm_inplace, HcLikelihoodEngineConfig};
use rust_htslib::bam;
use rust_htslib::bam::record::{Aux, Cigar, CigarString};

const MAPQ: u8 = 25;
const MIN_TAIL: u8 = 9;
const RAW_SEQ: &[u8] = b"AAAAACGT";
const RAW_BQ: &[u8] = &[31, 5, 40, 31, 31, 31, 31, 2];
const RAW_BI: &[u8] = &[45, 45, 45, 45, 30, 22, 30, 33];
const RAW_BD: &[u8] = &[40, 40, 40, 40, 25, 22, 25, 28];

#[derive(Clone, Debug, Eq, PartialEq)]
struct KernelTuple {
    read_id: String,
    hap_id: String,
    read_bases: Vec<u8>,
    hap_bases: Vec<u8>,
    bq: Vec<u8>,
    iq: Vec<u8>,
    dq: Vec<u8>,
    gcp: Vec<u8>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum ArrayName {
    ReadBases,
    HapBases,
    Bq,
    Iq,
    Dq,
    Gcp,
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct ElemDiff {
    array: ArrayName,
    index: usize,
    java: u8,
    rust: u8,
}

fn first_elem_diff(java: &KernelTuple, rust: &KernelTuple) -> Option<ElemDiff> {
    for (array, ja, ra) in [
        (
            ArrayName::ReadBases,
            java.read_bases.as_slice(),
            rust.read_bases.as_slice(),
        ),
        (
            ArrayName::HapBases,
            java.hap_bases.as_slice(),
            rust.hap_bases.as_slice(),
        ),
        (ArrayName::Bq, java.bq.as_slice(), rust.bq.as_slice()),
        (ArrayName::Iq, java.iq.as_slice(), rust.iq.as_slice()),
        (ArrayName::Dq, java.dq.as_slice(), rust.dq.as_slice()),
        (ArrayName::Gcp, java.gcp.as_slice(), rust.gcp.as_slice()),
    ] {
        if ja.len() != ra.len() {
            return Some(ElemDiff {
                array,
                index: ja.len().min(ra.len()),
                java: ja.len() as u8,
                rust: ra.len() as u8,
            });
        }
        if let Some(i) = ja.iter().zip(ra.iter()).position(|(a, b)| a != b) {
            return Some(ElemDiff {
                array,
                index: i,
                java: ja[i],
                rust: ra[i],
            });
        }
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

fn make_read(qname: &[u8], seq: &[u8], bq: &[u8], bi: &[u8], bd: &[u8]) -> bam::Record {
    let mut rec = bam::Record::new();
    let cigar = CigarString::from(vec![Cigar::Match(seq.len() as u32)]);
    rec.set(qname, Some(&cigar), seq, bq);
    rec.set_mapq(MAPQ);
    rec.set_pos(100);
    rec.set_tid(0);
    rec.push_aux(b"BI", Aux::String(&phred_to_fastq(bi)))
        .expect("BI");
    rec.push_aux(b"BD", Aux::String(&phred_to_fastq(bd)))
        .expect("BD");
    rec
}

/// Java values immediately before `subComputeReadLikelihoodGivenHaplotypeLog10`.
fn java_kernel_tuple(
    read_id: &str,
    hap_id: &str,
    read_bases: &[u8],
    raw_bq: &[u8],
    mapq: u8,
    bi: Option<&[u8]>,
    bd: Option<&[u8]>,
    hap_bases: &[u8],
) -> KernelTuple {
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
    KernelTuple {
        read_id: read_id.to_string(),
        hap_id: hap_id.to_string(),
        read_bases: read_bases.to_vec(),
        hap_bases: hap_bases.to_vec(),
        bq,
        iq,
        dq,
        gcp: vec![GATK_PARITY_DEFAULT_GCP; read_bases.len()],
    }
}

/// Production `score_read_against_haplotypes` planes (no IQ/DQ floor).
fn rust_kernel_tuple_from_record(
    read_id: &str,
    hap_id: &str,
    rec: &bam::Record,
    hap_bases: &[u8],
) -> KernelTuple {
    let cfg = HcLikelihoodEngineConfig::default();
    let read_bases = rec.seq().as_bytes();
    let mut bq = rec.qual().to_vec();
    prepare_read_quals_for_pairhmm_inplace(&mut bq, rec.mapq(), &cfg);
    let bi = aux_phred(rec, b"BI");
    let bd = aux_phred(rec, b"BD");
    let mut iq = indel_gop_from_optional_tag(bi.as_deref(), read_bases.len()).unwrap();
    let mut dq = indel_gop_from_optional_tag(bd.as_deref(), read_bases.len()).unwrap();
    apply_pcr_error_model(&read_bases, &mut iq, &mut dq, PcrErrorModel::Conservative);
    KernelTuple {
        read_id: read_id.to_string(),
        hap_id: hap_id.to_string(),
        read_bases,
        hap_bases: hap_bases.to_vec(),
        bq,
        iq,
        dq,
        gcp: vec![GATK_PARITY_DEFAULT_GCP; rec.seq_len()],
    }
}

fn rust_eligible(haps: &[Haplotype]) -> Vec<usize> {
    let ref_len = haps
        .iter()
        .find(|h| h.is_reference)
        .map(|h| h.bases.len())
        .unwrap_or(0);
    if ref_len == 0 {
        return (0..haps.len()).collect();
    }
    let max_alt = ref_len.saturating_add(8);
    haps.iter()
        .enumerate()
        .filter(|(_, h)| h.is_reference || h.bases.len() <= max_alt)
        .map(|(i, _)| i)
        .collect()
}

fn clipped_fixture_read(qname: &[u8]) -> bam::Record {
    let rec = make_read(qname, RAW_SEQ, RAW_BQ, RAW_BI, RAW_BD);
    hard_clip_low_qual_ends(&rec, MIN_TAIL)
}

fn fixture_haps() -> Vec<Haplotype> {
    let clipped = clipped_fixture_read(b"r0");
    let bases = clipped.seq().as_bytes();
    let mut alt = bases.clone();
    alt[3] = b'T';
    let mut haps = vec![Haplotype::new(alt, false), Haplotype::new(bases, true)];
    haps[0].score = 10.0;
    haps[1].score = 1.0;
    sort_haplotypes_assembly_result_order(&mut haps);
    haps
}

fn pair_matrix(
    reads: &[bam::Record],
    haps: &[Haplotype],
) -> Vec<(usize, usize, KernelTuple, KernelTuple)> {
    let eligible = rust_eligible(haps);
    let mut out = Vec::new();
    for (ri, rec) in reads.iter().enumerate() {
        let qname = String::from_utf8_lossy(rec.qname()).into_owned();
        for &hi in &eligible {
            let hap = &haps[hi];
            let hap_id = format!(
                "{}:{}",
                if hap.is_reference { "REF" } else { "ALT" },
                String::from_utf8_lossy(&hap.bases)
            );
            let rust = rust_kernel_tuple_from_record(&qname, &hap_id, rec, &hap.bases);
            let java = java_kernel_tuple(
                &qname,
                &hap_id,
                &rec.seq().as_bytes(),
                rec.qual(),
                rec.mapq(),
                aux_phred(rec, b"BI").as_deref(),
                aux_phred(rec, b"BD").as_deref(),
                &hap.bases,
            );
            out.push((ri, hi, java, rust));
        }
    }
    out
}

#[test]
fn complete_tuple_is_elementwise_identical_on_clipped_bi_bd_pcr_fixture() {
    let reads = vec![clipped_fixture_read(b"r0"), clipped_fixture_read(b"r1")];
    assert_eq!(reads[0].seq().as_bytes(), b"AAAAACG");
    let haps = fixture_haps();
    assert_eq!(haps.len(), 2);
    assert!(!haps[0].is_reference);
    assert!(haps[1].is_reference);

    let pairs = pair_matrix(&reads, &haps);
    assert_eq!(pairs.len(), 4, "2 reads × 2 haplotypes");
    assert_eq!(
        pairs
            .iter()
            .map(|(r, h, _, _)| (*r, *h))
            .collect::<Vec<_>>(),
        vec![(0, 0), (0, 1), (1, 0), (1, 1)]
    );

    eprintln!(
        "6R.75 fixture first pair: read={} hap={} seq={} bq={:?} iq={:?} dq={:?} gcp={:?}",
        pairs[0].2.read_id,
        pairs[0].2.hap_id,
        String::from_utf8_lossy(&pairs[0].2.read_bases),
        pairs[0].2.bq,
        pairs[0].2.iq,
        pairs[0].2.dq,
        pairs[0].2.gcp
    );
    for (ri, hi, java, rust) in &pairs {
        assert_eq!(
            first_elem_diff(java, rust),
            None,
            "pair ({ri},{hi}) first diff java_id={} rust_id={} java={java:?} rust={rust:?}",
            java.read_id,
            rust.read_id
        );
        assert_eq!(java.read_bases.len(), java.bq.len());
        assert_eq!(java.bq.len(), java.iq.len());
        assert_eq!(java.iq.len(), java.dq.len());
        assert_eq!(java.dq.len(), java.gcp.len());
        assert_eq!(java.hap_bases.len(), rust.hap_bases.len());
    }
}

#[test]
fn quality_arrays_are_nontrivial_and_not_pcr_masked_to_q45() {
    let rec = clipped_fixture_read(b"r0");
    let hap = rec.seq().as_bytes();
    let rust = rust_kernel_tuple_from_record("r0", "REF", &rec, &hap);
    assert_eq!(rust.bq[1], 6, "BQ 5 is capped to MIN_USABLE");
    assert_eq!(rust.bq[2], 25, "BQ 40 is capped by MAPQ 25");
    assert!(
        rust.iq.iter().any(|&q| q != 45 && q != 40),
        "IQ must not be a flat Q45/Q40 mask: {:?}",
        rust.iq
    );
    assert!(
        rust.dq.iter().any(|&q| q != 45 && q != 40),
        "DQ must not be a flat Q45/Q40 mask: {:?}",
        rust.dq
    );
    assert!(rust.gcp.iter().all(|&q| q == 10));
    assert!(
        rust.iq
            .iter()
            .any(|&q| q != rust.dq.iter().copied().next().unwrap_or(0))
            || rust.iq != rust.dq,
        "IQ and DQ are independent: iq={:?} dq={:?}",
        rust.iq,
        rust.dq
    );
    assert!(
        rust.iq.iter().all(|&q| q >= MIN_USABLE_Q_SCORE)
            && rust.dq.iter().all(|&q| q >= MIN_USABLE_Q_SCORE),
        "IQ/DQ floor is a no-op on this fixture"
    );
}

#[test]
fn pcr_adjusts_at_least_one_iq_cell_below_raw_bi() {
    let rec = clipped_fixture_read(b"r0");
    let bi = aux_phred(&rec, b"BI").expect("BI");
    let hap = rec.seq().as_bytes();
    let rust = rust_kernel_tuple_from_record("r0", "REF", &rec, &hap);
    let lowered = rust
        .iq
        .iter()
        .zip(bi.iter())
        .filter(|(iq, raw)| *iq < *raw)
        .count();
    assert!(
        lowered > 0,
        "PCR must write at least one IQ cell: iq={:?} bi={:?}",
        rust.iq,
        bi
    );
    assert!(
        rust.iq.iter().zip(bi.iter()).any(|(iq, raw)| iq == raw),
        "at least one IQ cell must remain raw BI (not Q45): iq={:?} bi={:?}",
        rust.iq,
        bi
    );
}

/// Anti-mask: Java `callRegion` drops `unclippedReadLength < 10` after trim.
/// Before 6R.75, `clip_finalized_reads_in_place` kept 2/5/7/9 bp remnants and PairHMM scored them.
#[test]
fn post_trim_stubs_shorter_than_ten_are_dropped_before_pairhmm() {
    use gatk_haplotypecaller::assembly_region_finalize::clip_finalized_reads_to_region;
    use gatk_haplotypecaller::assembly_region_iterator::AssemblyRegion;
    use gatk_haplotypecaller::{FeatureContext, GenomePosition, ReferenceContext};

    fn region() -> AssemblyRegion {
        AssemblyRegion {
            contig: "20".into(),
            start: GenomePosition::new_1based(101),
            end: GenomePosition::new_1based(200),
            extended_start: GenomePosition::new_1based(101),
            extended_end: GenomePosition::new_1based(200),
            extension: 0,
            reads: vec![],
            read_qnames: vec![],
            reference: ReferenceContext::empty(),
            features: FeatureContext::empty(),
            pileup_loci: vec![],
            is_active: true,
        }
    }
    let mut stub = make_read(b"stub7", b"CATGGAG", &[30; 7], &[45; 7], &[45; 7]);
    stub.set_flags(0);
    assert_eq!(stub.seq_len(), 7);
    let kept_stub = clip_finalized_reads_to_region(&[stub], &region());
    assert!(
        kept_stub.is_empty(),
        "anti-mask: 7 bp remnant must not reach PairHMM: {:?}",
        kept_stub.iter().map(|r| r.seq_len()).collect::<Vec<_>>()
    );
    let mut ok = make_read(b"ok10", b"CATGGAGCCG", &[30; 10], &[45; 10], &[45; 10]);
    ok.set_flags(0);
    let kept_ok = clip_finalized_reads_to_region(&[ok], &region());
    assert_eq!(kept_ok.len(), 1);
    assert_eq!(kept_ok[0].seq_len(), 10);
}

#[test]
fn haplotype_selection_is_full_assembly_list_in_java_order() {
    let haps = fixture_haps();
    let eligible = rust_eligible(&haps);
    assert_eq!(eligible, vec![0, 1]);
    assert_eq!(haps[0].bases, b"AAATACG");
    assert_eq!(haps[1].bases, b"AAAAACG");
    let ids: Vec<String> = haps
        .iter()
        .map(|h| format!("{}:{}", h.is_reference, String::from_utf8_lossy(&h.bases)))
        .collect();
    assert_eq!(
        ids,
        vec!["false:AAATACG".to_string(), "true:AAAAACG".to_string()]
    );
}

#[test]
fn read_selection_uses_qname_order_of_the_clipped_list() {
    let reads = vec![clipped_fixture_read(b"r0"), clipped_fixture_read(b"r1")];
    let names: Vec<String> = reads
        .iter()
        .map(|r| String::from_utf8_lossy(r.qname()).into_owned())
        .collect();
    assert_eq!(names, vec!["r0".to_string(), "r1".to_string()]);
    assert_eq!(reads.len(), 2);
}

#[test]
fn kernel_receives_unpadded_sequences_gcp_length_equals_read() {
    let rec = clipped_fixture_read(b"r0");
    let hap = rec.seq().as_bytes();
    let rust = rust_kernel_tuple_from_record("r0", "REF", &rec, &hap);
    assert_eq!(rust.read_bases.len(), 7);
    assert_eq!(rust.hap_bases.len(), 7);
    assert_eq!(rust.gcp.len(), 7);
    assert!(rust
        .read_bases
        .iter()
        .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'N')));
}

#[derive(Clone, Debug)]
struct DumpRegion {
    n_reads: usize,
    n_haps: usize,
    haps: Vec<String>,
    /// Processed-read identity in dump order: (bases, bq, iq, dq, gcp).
    reads: Vec<(String, String, String, String, String)>,
    rows: Vec<(String, String, String, String, String, String)>,
}

fn parse_pairhmm_dump(text: &str) -> Vec<DumpRegion> {
    let mut rows = Vec::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split(' ').collect();
        if parts.len() < 6 {
            continue;
        }
        rows.push((
            parts[0].to_string(),
            parts[1].to_string(),
            parts[2].to_string(),
            parts[3].to_string(),
            parts[4].to_string(),
            parts[5].to_string(),
        ));
    }
    let mut regions = Vec::new();
    let mut i = 0;
    while i < rows.len() {
        let rid = (
            rows[i].1.clone(),
            rows[i].2.clone(),
            rows[i].3.clone(),
            rows[i].4.clone(),
            rows[i].5.clone(),
        );
        let mut haps = Vec::new();
        let mut j = i;
        while j < rows.len() {
            let jr = (
                rows[j].1.clone(),
                rows[j].2.clone(),
                rows[j].3.clone(),
                rows[j].4.clone(),
                rows[j].5.clone(),
            );
            if jr != rid {
                break;
            }
            haps.push(rows[j].0.clone());
            j += 1;
        }
        let h = haps.len();
        let mut reads = vec![rid];
        let mut k = j;
        while k + h <= rows.len() {
            let block_haps: Vec<String> = (0..h).map(|t| rows[k + t].0.clone()).collect();
            if block_haps != haps {
                break;
            }
            reads.push((
                rows[k].1.clone(),
                rows[k].2.clone(),
                rows[k].3.clone(),
                rows[k].4.clone(),
                rows[k].5.clone(),
            ));
            k += h;
        }
        regions.push(DumpRegion {
            n_reads: reads.len(),
            n_haps: h,
            haps,
            reads,
            rows: rows[i..k].to_vec(),
        });
        i = k;
    }
    regions
}

fn snp_motif_region(regions: &[DumpRegion]) -> Option<&DumpRegion> {
    regions
        .iter()
        .find(|r| r.haps.iter().any(|h| h.contains("GTGGCTCACGTCTGTAAT")))
}

fn first_fastq_diff(java: &str, rust: &str, array: ArrayName) -> Option<ElemDiff> {
    let j: Vec<u8> = java.bytes().map(|b| b.saturating_sub(33)).collect();
    let r: Vec<u8> = rust.bytes().map(|b| b.saturating_sub(33)).collect();
    if j.len() != r.len() {
        return Some(ElemDiff {
            array,
            index: j.len().min(r.len()),
            java: j.len() as u8,
            rust: r.len() as u8,
        });
    }
    j.iter()
        .zip(r.iter())
        .position(|(a, b)| a != b)
        .map(|i| ElemDiff {
            array,
            index: i,
            java: j[i],
            rust: r[i],
        })
}

/// Live capture for the report. Not the sole proof.
///
/// ```text
/// HOLDOUT_6R75=1 cargo test -p gatk-haplotypecaller --test forensic_6r75_pairhmm_complete_input_contract canonical_live -- --nocapture
/// ```
#[test]
fn canonical_live_complete_tuple_and_dimensions() {
    if std::env::var("HOLDOUT_6R75").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R75=1");
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
    use std::collections::{HashMap, HashSet};
    use std::path::Path;

    const INTERVAL: &str = "20:29455000-29456500";
    const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
    const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
    const JAVA_DUMP_REL: &str = "parity/giab/runs/local-pairhmm-diff/6r75_java_pairhmm_inputs.txt";
    const POS_SNP: u64 = 29_456_344;

    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    let java_dump_path = root.join(JAVA_DUMP_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: live BAM/ref missing");
        return;
    }
    let rust_dump_path = std::env::temp_dir().join("6r75_rust_pairhmm_inputs.txt");
    std::env::set_var(
        "GATK_RS_PAIRHMM_INPUT_DUMP",
        rust_dump_path.to_string_lossy().as_ref(),
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
    let eligible = rust_eligible(&outcome.assembly.haplotypes);
    assert_eq!(
        eligible.len(),
        outcome.assembly.haplotypes.len(),
        "ref+8 eligibility must be a no-op on this assembly"
    );
    assert_eq!(
        eligible,
        (0..outcome.assembly.haplotypes.len()).collect::<Vec<_>>()
    );

    let rust_dump = std::fs::read_to_string(&rust_dump_path).expect("rust dump");
    let rust_regions = parse_pairhmm_dump(&rust_dump);
    let rust_reg = snp_motif_region(&rust_regions).unwrap_or_else(|| {
        panic!(
            "no SNP-motif region in rust dump; n_regions={} sizes={:?}",
            rust_regions.len(),
            rust_regions
                .iter()
                .map(|r| (r.n_reads, r.n_haps, r.haps.first().map(|h| h.len())))
                .collect::<Vec<_>>()
        )
    });

    let mut first_java_cmp: Option<String> = None;
    if java_dump_path.is_file() {
        let java_dump = std::fs::read_to_string(&java_dump_path).expect("java dump");
        let java_regions = parse_pairhmm_dump(&java_dump);
        let java_reg = snp_motif_region(&java_regions).expect("java SNP-motif region");
        eprintln!(
            "6R.75 kernel-boundary: java {}×{} rust_dump {}×{} postfilter_ll {}×{} assembly_haps={}",
            java_reg.n_reads,
            java_reg.n_haps,
            rust_reg.n_reads,
            rust_reg.n_haps,
            n_reads.len(),
            n_haps.len(),
            outcome.assembly.haplotypes.len()
        );
        if java_reg.n_reads != rust_reg.n_reads {
            let java_set: HashSet<&String> = java_reg.reads.iter().map(|r| &r.0).collect();
            let rust_set: HashSet<&String> = rust_reg.reads.iter().map(|r| &r.0).collect();
            first_java_cmp = Some(format!(
                "read selection/order: java_n={} rust_n={} only_java_seq={} only_rust_seq={}",
                java_reg.n_reads,
                rust_reg.n_reads,
                java_set.difference(&rust_set).count(),
                rust_set.difference(&java_set).count()
            ));
        } else {
            let java_set: HashSet<&String> = java_reg.haps.iter().collect();
            let rust_set: HashSet<&String> = rust_reg.haps.iter().collect();
            let read_order_diff = java_reg.reads.iter().map(|r| &r.0).collect::<Vec<_>>()
                != rust_reg.reads.iter().map(|r| &r.0).collect::<Vec<_>>();
            if java_reg.n_haps != rust_reg.n_haps {
                first_java_cmp = Some(format!(
                    "haplotype selection/order: java_n={} rust_n={} only_java={} only_rust={}",
                    java_reg.n_haps,
                    rust_reg.n_haps,
                    java_set.difference(&rust_set).count(),
                    rust_set.difference(&java_set).count()
                ));
            } else if java_set != rust_set {
                first_java_cmp = Some(format!(
                    "haplotype sequence: only_java={} only_rust={}",
                    java_set.difference(&rust_set).count(),
                    rust_set.difference(&java_set).count()
                ));
            } else if java_reg.haps != rust_reg.haps {
                first_java_cmp =
                    Some("haplotype selection/order: same set, different order".into());
            } else if read_order_diff {
                first_java_cmp =
                    Some("read selection/order: counts equal, sequence order differs".into());
            } else {
                let mut by_pair: HashMap<(String, String), (String, String, String, String)> =
                    HashMap::new();
                for (hap, read, bq, iq, dq, gcp) in &java_reg.rows {
                    by_pair.insert(
                        (read.clone(), hap.clone()),
                        (bq.clone(), iq.clone(), dq.clone(), gcp.clone()),
                    );
                }
                'pairs: for (hap, read, bq, iq, dq, gcp) in &rust_reg.rows {
                    let Some(j) = by_pair.get(&(read.clone(), hap.clone())) else {
                        first_java_cmp =
                            Some("read/haplotype sequence: rust pair missing in java".into());
                        break;
                    };
                    for (array, ja, ra) in [
                        (ArrayName::Bq, j.0.as_str(), bq.as_str()),
                        (ArrayName::Iq, j.1.as_str(), iq.as_str()),
                        (ArrayName::Dq, j.2.as_str(), dq.as_str()),
                        (ArrayName::Gcp, j.3.as_str(), gcp.as_str()),
                    ] {
                        if let Some(d) = first_fastq_diff(ja, ra, array) {
                            first_java_cmp = Some(format!(
                                "{:?} array, index {} java={} rust={}",
                                d.array, d.index, d.java, d.rust
                            ));
                            break 'pairs;
                        }
                    }
                }
            }
            eprintln!(
                "6R.75 hap_set_equal={} hap_order_equal={} read_order_equal={}",
                java_set == rust_set,
                java_reg.haps == rust_reg.haps,
                !read_order_diff
            );
        }
        eprintln!("6R.75 first_java_vs_rust_kernel_diff={first_java_cmp:?}");
    } else {
        eprintln!(
            "skip java dump compare: {} missing",
            java_dump_path.display()
        );
    }

    let finalized = finalize_region_reads_for_assembly(
        &covering[0].reads,
        covering[0],
        true,
        gatk_min_tail_quality_for_assembly(10),
        false,
    );
    let dummy_hap = outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .map(|h| h.bases.clone())
        .unwrap_or_else(|| b"N".to_vec());
    let mut n_floor = 0usize;
    let mut first_tuple_diff: Option<(usize, ElemDiff)> = None;
    let mut n_iq_ne_q40 = 0usize;
    for (i, rec) in finalized.iter().enumerate() {
        let rust = rust_kernel_tuple_from_record(
            &String::from_utf8_lossy(rec.qname()),
            "REF",
            rec,
            &dummy_hap,
        );
        let java = java_kernel_tuple(
            &rust.read_id,
            "REF",
            &rust.read_bases,
            rec.qual(),
            rec.mapq(),
            aux_phred(rec, b"BI").as_deref(),
            aux_phred(rec, b"BD").as_deref(),
            &dummy_hap,
        );
        n_floor += rust
            .iq
            .iter()
            .chain(rust.dq.iter())
            .filter(|&&q| q < MIN_USABLE_Q_SCORE)
            .count();
        n_iq_ne_q40 += rust.iq.iter().filter(|&&q| q != 40 && q != 45).count();
        if first_tuple_diff.is_none() {
            if let Some(d) = first_elem_diff(&java, &rust) {
                first_tuple_diff = Some((i, d));
            }
        }
    }
    let live = take_colocated_merge_numerics();
    let numerics = live.iter().find(|n| n.loc == POS_SNP);
    eprintln!(
        "6R.75 live: PL={:?} AD={:?} QUAL={:?} geno_reads={} ll_reads={} ll_haps={} assembly_haps={} eligible={} cells={} overlap={:?} walker_finalized={} floor_lt6={} iq_not_q40_or_q45={} contract_tuple_diff={:?} kernel={}x{}",
        vcf.samples[0].pl,
        vcf.samples[0].ad,
        vcf.quality,
        outcome.genotyping_reads.len(),
        n_reads.len(),
        n_haps.len(),
        outcome.assembly.haplotypes.len(),
        eligible.len(),
        outcome.read_likelihoods.len(),
        numerics.map(|n| (n.n_reads, n.pool_sizes.clone())),
        finalized.len(),
        n_floor,
        n_iq_ne_q40,
        first_tuple_diff,
        rust_reg.n_reads,
        rust_reg.n_haps,
    );
    assert_eq!(n_reads.len(), 136);
    assert_eq!(n_haps.len(), 68);
    assert_eq!(first_tuple_diff, None, "prep-path Java vs Rust tuples");
    assert_eq!(n_floor, 0, "IQ/DQ floor is a no-op on this BAM");
    assert!(
        rust_reg.reads.iter().all(|r| r.0.len() >= 10),
        "PairHMM must not score Java post-trim stubs"
    );
    assert_eq!(
        rust_reg.n_reads, 153,
        "Java kernel read count after stub drop"
    );
    assert_eq!(
        rust_reg.n_haps, 70,
        "6R.81: production kernel haplotype count after trimDown ref+alt collapse"
    );
    if let Some(diff) = first_java_cmp {
        if diff.contains("java_n=")
            && diff.contains("rust_n=")
            && diff.starts_with("read selection")
        {
            panic!("FIRST REMAINING DIVERGENCE: {diff}");
        }
        eprintln!("6R.75 next remaining after stub filter (not stacked): {diff}");
    }
}
