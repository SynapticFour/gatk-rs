//! 6R.231: first production operation that clips the five poorly-modeled
//! reads at `20:29455379 G/A` from the Java padded-window representation
//! to start `29455355`.
//!
//! Proof-only. Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE.
//!
//! 6R.230 closed `READ_SEQUENCE_INPUT_DIVERGENCE`. This round traces the
//! clip/trim lifecycle. Haplotype population is deferred unless clipping
//! depends on a haplotype-derived interval.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r231_read_clipping_boundary -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly_based_caller::assemble_reads_with_finalized;
use gatk_haplotypecaller::assembly_region_finalize::{
    clip_finalized_reads_to_region, finalize_region_reads_for_assembly,
    gatk_min_tail_quality_for_assembly,
};
use gatk_haplotypecaller::assembly_region_trimmer::{AssemblyRegionTrimmer, TrimVariant};
use gatk_haplotypecaller::read_unclip::{alignment_end_1based, hard_clip_to_region};
use gatk_haplotypecaller::reference_context::ReferenceContext;
use gatk_haplotypecaller::{
    begin_hap_list_observe, begin_realign_observe, call_disposition, flatten_assembly_regions,
    take_hap_list_trim_span, take_realign_observe, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams,
    WalkerTraversalConfig,
};
use rust_htslib::bam::record::Cigar;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_455_379;
const JAVA_VARIANT_START: u64 = 29_455_314;
const JAVA_VARIANT_END: u64 = 29_455_523;
const JAVA_PADDED_START: u64 = 29_455_294;
const JAVA_PADDED_END: u64 = 29_455_584;
const SNP_PADDING: u64 = 20;
const RUST_CLIP_START: u64 = 29_455_294;
const FIVE: &[(&str, u16)] = &[
    ("HISEQ1:13:H8G92ADXX:1:2102:7192:18079", 99),
    ("HISEQ1:9:H8962ADXX:1:1116:1789:43193", 163),
    ("HWI-D00360:6:H81VLADXX:1:1111:8050:25694", 163),
    ("HWI-D00360:6:H81VLADXX:1:2102:1733:39463", 163),
    ("HWI-D00360:7:H88WKADXX:2:1201:4043:96748", 163),
];
const JAVA_PAIRHMM: &[(&str, u16, i64, usize, &str)] = &[
    (
        "HISEQ1:13:H8G92ADXX:1:2102:7192:18079",
        99,
        29_455_326,
        148,
        "148M",
    ),
    (
        "HISEQ1:9:H8962ADXX:1:1116:1789:43193",
        163,
        29_455_294,
        113,
        "33H113M2H",
    ),
    (
        "HWI-D00360:6:H81VLADXX:1:1111:8050:25694",
        163,
        29_455_328,
        148,
        "148M",
    ),
    (
        "HWI-D00360:6:H81VLADXX:1:2102:1733:39463",
        163,
        29_455_319,
        148,
        "148M",
    ),
    (
        "HWI-D00360:7:H88WKADXX:2:1201:4043:96748",
        163,
        29_455_312,
        148,
        "148M",
    ),
];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R231\t{key}\t{}", value.as_ref());
}

fn fnv1a64_hex(data: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in data {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{h:016x}")
}

fn cigar_str(rec: &rust_htslib::bam::Record) -> String {
    format!("{}", rec.cigar())
}

fn leading_hard(rec: &rust_htslib::bam::Record) -> u32 {
    rec.cigar()
        .iter()
        .take_while(|c| matches!(c, Cigar::HardClip(_)))
        .map(|c| match c {
            Cigar::HardClip(n) => *n,
            _ => 0,
        })
        .sum()
}

fn leading_soft(rec: &rust_htslib::bam::Record) -> u32 {
    rec.cigar()
        .iter()
        .skip_while(|c| matches!(c, Cigar::HardClip(_)))
        .take_while(|c| matches!(c, Cigar::SoftClip(_)))
        .map(|c| match c {
            Cigar::SoftClip(n) => *n,
            _ => 0,
        })
        .sum()
}

fn trailing_hard(rec: &rust_htslib::bam::Record) -> u32 {
    rec.cigar()
        .iter()
        .rev()
        .take_while(|c| matches!(c, Cigar::HardClip(_)))
        .map(|c| match c {
            Cigar::HardClip(n) => *n,
            _ => 0,
        })
        .sum()
}

fn trailing_soft(rec: &rust_htslib::bam::Record) -> u32 {
    rec.cigar()
        .iter()
        .rev()
        .skip_while(|c| matches!(c, Cigar::HardClip(_)))
        .take_while(|c| matches!(c, Cigar::SoftClip(_)))
        .map(|c| match c {
            Cigar::SoftClip(n) => *n,
            _ => 0,
        })
        .sum()
}

#[derive(Clone, Debug)]
struct Snap {
    stage: &'static str,
    start: i64,
    end: i64,
    len: usize,
    cigar: String,
    seq_hash: String,
    bq_hash: String,
    seq: Vec<u8>,
    bq: Vec<u8>,
    hard_left: u32,
    hard_right: u32,
    soft_left: u32,
    soft_right: u32,
}

fn snap_rec(stage: &'static str, rec: &rust_htslib::bam::Record) -> Snap {
    let seq = rec.seq().as_bytes();
    let bq = rec.qual().to_vec();
    Snap {
        stage,
        start: rec.pos() + 1,
        end: i64::from(alignment_end_1based(rec)),
        len: seq.len(),
        cigar: cigar_str(rec),
        seq_hash: fnv1a64_hex(&seq),
        bq_hash: fnv1a64_hex(&bq),
        seq,
        bq,
        hard_left: leading_hard(rec),
        hard_right: trailing_hard(rec),
        soft_left: leading_soft(rec),
        soft_right: trailing_soft(rec),
    }
}

fn dump_snap(q: &str, f: u16, s: &Snap, orig_start: i64, orig_end: i64, seq: bool) {
    kv(
        "snap",
        format!(
            "qname={q}\tFLAG={f}\tcontig=20\tstage={}\tstart={}\tend={}\tseq_len={}\tcigar={}\thard_clip_left={}\thard_clip_right={}\tsoft_clip_left={}\tsoft_clip_right={}\tseq_hash={}\tbq_hash={}\toriginal_start={orig_start}\toriginal_end={orig_end}",
            s.stage,
            s.start,
            s.end,
            s.len,
            s.cigar,
            s.hard_left,
            s.hard_right,
            s.soft_left,
            s.soft_right,
            s.seq_hash,
            s.bq_hash
        ),
    );
    if seq {
        kv(
            "snap_seq",
            format!(
                "qname={q}\tstage={}\tseq={}",
                s.stage,
                String::from_utf8_lossy(&s.seq)
            ),
        );
    }
}

fn find_five<'a>(
    reads: impl IntoIterator<Item = &'a rust_htslib::bam::Record>,
    q: &str,
    f: u16,
) -> Option<&'a rust_htslib::bam::Record> {
    reads
        .into_iter()
        .find(|r| String::from_utf8_lossy(r.qname()) == q && r.flags() == f)
}

#[test]
fn forensic_6r231_read_clipping_boundary() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv("target", "20:29455379 G/A");
    kv(
        "predecessor",
        "6R.230 READ_SEQUENCE_INPUT_DIVERGENCE CLOSED",
    );
    kv(
        "preserve",
        "6R.226 cache; 6R.227 genotype_from_marginalized_rows; 6R.228 retainEvidence; 6R.229 poorly-modeled predicate; 6R.230 PairHMM inputs unpatched",
    );

    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        panic!("missing chr20_tiny BAM/REF");
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
    let region = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("ActiveFull");

    kv(
        "active_region",
        format!("20:{}-{}", region.start.get(), region.end.get()),
    );
    kv(
        "assembly_region_padded",
        format!(
            "20:{}-{}",
            region.extended_start.get(),
            region.extended_end.get()
        ),
    );
    kv(
        "java_variant_interval",
        format!("20:{JAVA_VARIANT_START}-{JAVA_VARIANT_END}"),
    );
    kv(
        "java_padded_variant_interval",
        format!("20:{JAVA_PADDED_START}-{JAVA_PADDED_END}"),
    );
    kv(
        "java_pairhmm_read_interval",
        format!("20:{JAVA_PADDED_START}-{JAVA_PADDED_END}"),
    );
    assert_eq!(JAVA_PADDED_START, JAVA_VARIANT_START - SNP_PADDING);
    kv(
        "java_padded_provenance",
        format!("{JAVA_PADDED_START} = {JAVA_VARIANT_START} - snp_padding({SNP_PADDING})"),
    );
    let java_stages = include_str!("forensic_6r229_java_stages.tsv");
    assert!(
        java_stages.contains("variant_span\t20:29455314-29455523"),
        "6R.229 Java dump pins variant_span 20:29455314-29455523"
    );
    assert!(
        java_stages.contains("variant_padded\t20:29455294-29455584"),
        "6R.229 Java dump pins variant_padded 20:29455294-29455584"
    );

    let args = CallRegionArgs::strict_java();
    let mut region_for_assemble = region.clone();
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let assembled = assemble_reads_with_finalized(
        &mut region_for_assemble,
        &dict,
        &mut ref_cache,
        &args.assemble,
    )
    .expect("assemble");
    let mut trim_variants: Vec<TrimVariant> = assembled
        .assembly
        .variation_events()
        .iter()
        .map(|e| TrimVariant {
            contig: e.contig.clone(),
            start: e.start_1based.get(),
            end: e.end_1based.get(),
            is_indel: e.is_indel(),
            ref_allele: e.ref_allele.clone(),
            alt_allele: e.alt_allele.clone(),
        })
        .collect();
    trim_variants.sort_by_key(|v| (v.start, v.end, v.is_indel));
    let overlapping: Vec<&TrimVariant> = trim_variants
        .iter()
        .filter(|v| v.overlaps_active_region(region))
        .collect();
    kv("trim_variant_n", trim_variants.len().to_string());
    kv("trim_variant_overlapping_n", overlapping.len().to_string());
    for v in &overlapping {
        kv(
            "trim_event",
            format!(
                "start={}\tend={}\tindel={}\tcontig={}",
                v.start, v.end, v.is_indel, v.contig
            ),
        );
    }
    let has_java_leftmost = trim_variants.iter().any(|v| v.start == JAVA_VARIANT_START);
    kv(
        "rust_variation_events_contain_java_leftmost_29455314",
        has_java_leftmost.to_string(),
    );
    assert!(
        has_java_leftmost,
        "6R.239: EventMap contains Java leftmost SNP 29455314 G>C"
    );
    let leftmost = overlapping
        .iter()
        .min_by_key(|v| v.start)
        .expect("at least one overlapping trim event");
    kv(
        "leftmost_trim_event",
        format!(
            "start={}\tend={}\tindel={}",
            leftmost.start, leftmost.end, leftmost.is_indel
        ),
    );
    assert_eq!(leftmost.start, JAVA_VARIANT_START);
    assert!(
        !leftmost.is_indel,
        "leftmost overlapping trim event is a SNP so padding is 20"
    );
    let expected_clip = leftmost.start.saturating_sub(SNP_PADDING);
    kv(
        "rust_clip_provenance",
        format!(
            "{RUST_CLIP_START} = leftmost_trim_snp({}) - snp_padding({SNP_PADDING}) = {expected_clip}",
            leftmost.start
        ),
    );
    assert_eq!(expected_clip, RUST_CLIP_START);

    let trimmer = AssemblyRegionTrimmer::new(args.trimmer.clone(), &dict, &region.contig);
    let mut ref_cache2 = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let ref_ctx = ReferenceContext::from_interval(
        &dict,
        &mut ref_cache2,
        &region.contig,
        region.extended_start.get(),
        region.extended_end.get(),
    )
    .expect("ref_ctx");
    let trim_result = trimmer.trim(region, &trim_variants, Some(&ref_ctx));
    assert!(trim_result.variation_present);
    kv(
        "trim_result_variant",
        format!(
            "{}-{}",
            trim_result.variant_start.unwrap(),
            trim_result.variant_end.unwrap()
        ),
    );
    kv(
        "trim_result_padded",
        format!(
            "{}-{}",
            trim_result.padded_variant_start.unwrap(),
            trim_result.padded_variant_end.unwrap()
        ),
    );
    assert_eq!(trim_result.padded_variant_start.unwrap(), RUST_CLIP_START);

    let region_for_gt = AssemblyRegionTrimmer::apply_trim(region, &trim_result);
    kv(
        "region_for_genotyping_active",
        format!(
            "20:{}-{}",
            region_for_gt.start.get(),
            region_for_gt.end.get()
        ),
    );
    kv(
        "region_for_genotyping_padded",
        format!(
            "20:{}-{}",
            region_for_gt.extended_start.get(),
            region_for_gt.extended_end.get()
        ),
    );
    assert_eq!(region_for_gt.extended_start.get(), RUST_CLIP_START);

    let finalized = finalize_region_reads_for_assembly(
        &region.reads,
        region,
        true,
        gatk_min_tail_quality_for_assembly(10),
        false,
    );
    let trim_clipped = clip_finalized_reads_to_region(&finalized, &region_for_gt);

    begin_hap_list_observe();
    begin_realign_observe();
    let _outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let span = take_hap_list_trim_span().expect("trim span");
    let realign = take_realign_observe();
    kv(
        "hap_list_trim_span",
        format!(
            "active={}-{} trim={}-{}",
            span.active_start, span.active_end, span.trim_start, span.trim_end
        ),
    );
    assert_eq!(span.trim_start, RUST_CLIP_START);
    kv(
        "first_interval_divergence",
        format!(
            "Java padded variant {JAVA_PADDED_START}-{JAVA_PADDED_END} vs Rust region_for_genotyping.extended {}-{}",
            region_for_gt.extended_start.get(),
            region_for_gt.extended_end.get()
        ),
    );
    kv(
        "java_controlling_interval",
        format!(
            "AssemblyRegionTrimmer padded variant span 20:{JAVA_PADDED_START}-{JAVA_PADDED_END}"
        ),
    );
    kv(
        "rust_controlling_interval",
        format!(
            "AssemblyRegionTrimmer::apply_trim → region.extended {}-{} then clip_finalized_reads_in_place/hard_clip_to_region",
            region_for_gt.extended_start.get(),
            region_for_gt.extended_end.get()
        ),
    );

    let mut first_ops = BTreeMap::new();
    for (i, (q, f)) in FIVE.iter().enumerate() {
        let orig = region
            .reads
            .iter()
            .find(|r| String::from_utf8_lossy(r.qname()) == *q && r.flags() == *f)
            .unwrap_or_else(|| panic!("orig missing {q}"));
        let after_fin =
            find_five(finalized.iter(), q, *f).unwrap_or_else(|| panic!("finalize missing {q}"));
        let after_trim_clip = find_five(trim_clipped.iter(), q, *f)
            .unwrap_or_else(|| panic!("trim-clip missing {q}"));
        let pairhmm = realign
            .iter()
            .find(|r| r.qname == *q && r.flags == *f)
            .unwrap_or_else(|| panic!("pairhmm snap missing {q}"));

        let s_orig = snap_rec("orig_bam", orig.as_ref());
        let s_fin = snap_rec("after_finalize_assembly_padded", after_fin);
        let s_trim = snap_rec("after_clip_finalized_reads_to_trim_span", after_trim_clip);
        assert_eq!(pairhmm.orig_start_1based, s_trim.start);
        assert_eq!(pairhmm.orig_cigar, s_trim.cigar);
        assert_eq!(pairhmm.orig_seq_len, s_trim.len);

        let java = JAVA_PAIRHMM[i];
        assert_eq!(java.0, *q);

        dump_snap(q, *f, &s_orig, s_orig.start, s_orig.end, false);
        dump_snap(q, *f, &s_fin, s_orig.start, s_orig.end, false);
        dump_snap(q, *f, &s_trim, s_orig.start, s_orig.end, true);

        let java_clip = hard_clip_to_region(after_fin, JAVA_PADDED_START, JAVA_PADDED_END);
        let s_java = snap_rec("java_hardClipToRegion_padded_variant", &java_clip);
        assert_eq!(s_java.start, java.2, "{q} Java counterpart clip start");
        assert_eq!(s_java.len, java.3, "{q} Java counterpart clip length");
        assert_eq!(s_java.cigar, java.4, "{q} Java counterpart clip CIGAR");
        dump_snap(q, *f, &s_java, s_orig.start, s_orig.end, true);

        kv(
            "lifecycle",
            format!(
                "qname={q}\tFLAG={f}\torig={}/{} {}\tfinalize={}/{} {}\tjava_clip={}/{} {}\ttrim_clip={}/{} {}\tpairhmm={}/{} {}\tjava_pin={}/{} {}",
                s_orig.start,
                s_orig.len,
                s_orig.cigar,
                s_fin.start,
                s_fin.len,
                s_fin.cigar,
                s_java.start,
                s_java.len,
                s_java.cigar,
                s_trim.start,
                s_trim.len,
                s_trim.cigar,
                pairhmm.orig_start_1based,
                pairhmm.orig_seq_len,
                pairhmm.orig_cigar,
                java.2,
                java.3,
                java.4
            ),
        );
        kv(
            "clip_case",
            format!(
                "qname={q}\tcase=A_physical_clip\tcoord_only=false\tsoft_left_orig={}\thard_left_fin={}\thard_left_java={}\thard_left_trim={}",
                s_orig.soft_left,
                s_fin.hard_left,
                s_java.hard_left,
                s_trim.hard_left
            ),
        );

        assert_ne!(
            s_fin.start, RUST_CLIP_START as i64,
            "{q} finalize to assembly padded span must not already start at {RUST_CLIP_START}"
        );
        assert_eq!(
            s_trim.start,
            java.2.max(RUST_CLIP_START as i64),
            "{q} trim clip start is Java-equivalent max(orig, padded_start)"
        );
        assert!(
            s_trim.len <= s_fin.len,
            "{q} trim clip does not lengthen the read"
        );
        let removed_from_fin = s_fin.len.saturating_sub(s_trim.len);
        if s_trim.len < s_fin.len {
            assert_eq!(
                &s_fin.seq[removed_from_fin..],
                s_trim.seq.as_slice(),
                "{q} Rust PairHMM sequence is an exact suffix of the post-finalize sequence"
            );
        }
        assert_eq!(
            &s_fin.bq[removed_from_fin..],
            s_trim.bq.as_slice(),
            "{q} BQ sliced by the same physical-clip offset"
        );
        assert!(
            s_java.len >= s_trim.len,
            "{q} Java padded-span clip is at least as long as Rust trim clip"
        );
        let removed_from_java = s_java.len - s_trim.len;
        assert_eq!(
            &s_java.seq[removed_from_java..],
            s_trim.seq.as_slice(),
            "{q} Rust PairHMM sequence is an exact suffix of the Java padded-window sequence"
        );
        assert_eq!(
            &s_java.bq[removed_from_java..],
            s_trim.bq.as_slice(),
            "{q} Rust BQ is the matching suffix of the Java padded-window BQ"
        );
        kv(
            "suffix",
            format!(
                "qname={q}\tprefix_removed_from_finalize={removed_from_fin}\tprefix_removed_from_java_clip={removed_from_java}\trust_is_suffix_of_finalize=true\trust_is_suffix_of_java=true\tinternal_diff=false\torientation=unchanged"
            ),
        );

        let last_java_eq = if s_fin.start == java.2 && s_fin.len == java.3 {
            "after_finalize_assembly_padded"
        } else {
            "after_finalize_assembly_padded (geometry still covers Java padded window left of 29455355); Java hardClipToRegion(29455294-29455584) is the Java-equivalent clip"
        };
        kv(
            "first_divergent_stage",
            format!(
                "qname={q}\tlast_java_equivalent={last_java_eq}\tfirst_divergent=clip_finalized_reads_in_place/hard_clip_to_region\tjava={}/{}\trust={}/{}",
                java.2, java.3, s_trim.start, s_trim.len
            ),
        );
        first_ops.insert(
            (*q).to_string(),
            "clip_finalized_reads_in_place / hard_clip_to_region(region.extended_start=29455355)",
        );
        assert_eq!(s_trim.start, java.2.max(RUST_CLIP_START as i64));
    }

    let shared = first_ops.values().next().copied().expect("ops");
    assert!(first_ops.values().all(|o| *o == shared));
    kv("shared_operation", shared);
    kv("shared_clipping_boundary", RUST_CLIP_START.to_string());
    kv(
        "haplotype_dependency",
        "trim interval is derived from untrimmed.variation_events (haplotype EventMap). Read clipping uses that interval via region.extended_start. Haplotype-population arrow is deferred; stop at this interval dependency.",
    );
    kv("classification", "READ_INTERVAL_INPUT_DIVERGENCE");
    kv(
        "first_causal_arrow",
        "AssemblyRegionTrimmer padded_variant_start = leftmost overlapping SNP − 20 = 29455355; clip_finalized_reads_in_place then physically hard-clips every frozen read to that interval. Java counterpart is ReadClipper.hardClipToRegion(paddedVariantSpan) with paddedVariantSpan=20:29455294-29455584.",
    );
    kv("pairhmm_downstream", "PairHMM consumes already-clipped evidence; 6R.230 likelihood split is downstream of this clip.");
    kv("production_src", "NONE");
}
