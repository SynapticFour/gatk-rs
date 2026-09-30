//! 6R.290: why one read is `51H97M` in Java and `60H88M` in Rust before PairHMM.
//! Localization only. No PL counterfactual.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r290_read_clipping_boundary -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::assembly_region_finalize::{
    finalize_region_reads_for_assembly, gatk_min_tail_quality_for_assembly,
};
use gatk_haplotypecaller::read_threading_assembler::DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH;
use gatk_haplotypecaller::read_unclip::{hard_clip_to_region, normalize_record_cigar};
use gatk_haplotypecaller::{
    begin_hap_list_observe, call_disposition, flatten_assembly_regions, take_hap_list_trim_span,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use rust_htslib::bam;
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_TSV_REL: &str = "gatk-haplotypecaller/tests/6r290_java_clip.tsv";
const QNAME: &str = "HISEQ1:11:H8GV6ADXX:1:2116:18670:99941";
const FLAGS: u16 = 99;
const TARGET: u64 = 29_455_649;
const EXTRA: &str = "CAAAGAGGA";

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R290\t{key}\t{}", value.as_ref());
}

fn find_read<'a>(reads: &'a [bam::Record], qname: &str, flags: u16) -> &'a bam::Record {
    reads
        .iter()
        .find(|r| r.flags() == flags && std::str::from_utf8(r.qname()).ok() == Some(qname))
        .expect("target read")
}

fn cigar_of(rec: &bam::Record) -> String {
    format!("{}", rec.cigar())
}

fn bases_of(rec: &bam::Record) -> Vec<u8> {
    rec.seq().as_bytes()
}

struct JavaRead {
    cigar: String,
    start: i64,
    end: i64,
    len: usize,
    bases: Vec<u8>,
}

fn kv_field(parts: &[String], key: &str) -> String {
    parts
        .iter()
        .find_map(|p| p.strip_prefix(&format!("{key}=")).map(|s| s.to_string()))
        .unwrap_or_default()
}

fn load_java(path: &Path) -> (String, String, JavaRead, JavaRead, JavaRead) {
    let text =
        std::fs::read_to_string(path).unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    let mut padded = String::new();
    let mut variant_padded = String::new();
    let mut iterator = None;
    let mut after_finalize = None;
    let mut after_trim = None;
    for line in text.lines() {
        let f: Vec<String> = line.split('\t').map(|s| s.to_string()).collect();
        if f.first().map(String::as_str) != Some("6R290") || f.len() < 3 {
            continue;
        }
        match f[1].as_str() {
            "padded_span" => padded = f[2].clone(),
            "variant_padded_span" => variant_padded = f[2].clone(),
            "iterator" => iterator = Some(parse_read(&f[2..])),
            "after_finalize" => after_finalize = Some(parse_read(&f[2..])),
            "after_trim_hardclip" => after_trim = Some(parse_read(&f[2..])),
            _ => {}
        }
    }
    (
        padded,
        variant_padded,
        iterator.expect("java iterator read"),
        after_finalize.expect("java after finalize"),
        after_trim.expect("java after trim hardclip"),
    )
}

fn parse_read(parts: &[String]) -> JavaRead {
    JavaRead {
        cigar: kv_field(parts, "cigar"),
        start: kv_field(parts, "start").parse().unwrap_or(0),
        end: kv_field(parts, "end").parse().unwrap_or(0),
        len: kv_field(parts, "len").parse().unwrap_or(0),
        bases: kv_field(parts, "bases").into_bytes(),
    }
}

#[test]
fn forensic_6r290_read_clipping_boundary() {
    assert_eq!(DEFAULT_NUM_BEST_HAPLOTYPES_PER_GRAPH, 128);
    let root = repo_root();
    let (java_padded, java_variant_padded, java_raw, java_final, java_trim) =
        load_java(&root.join(JAVA_TSV_REL));
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, INTERVAL).expect("interval");
    begin_hap_list_observe();
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
    let raw_owned: Vec<bam::Record> = region.reads.iter().map(|r| (**r).clone()).collect();
    let raw = find_read(&raw_owned, QNAME, FLAGS);
    let finalized = finalize_region_reads_for_assembly(
        &region.reads,
        region,
        true,
        gatk_min_tail_quality_for_assembly(10),
        false,
    );
    let after_finalize = find_read(&finalized, QNAME, FLAGS);
    let _outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let trim = take_hap_list_trim_span().expect("trim span");
    let mut clipped = hard_clip_to_region(after_finalize, trim.trim_start, trim.trim_end);
    normalize_record_cigar(&mut clipped);

    kv("qname", QNAME);
    kv("flags", FLAGS.to_string());
    kv(
        "rust_active",
        format!("{}-{}", region.start.get(), region.end.get()),
    );
    kv(
        "rust_padded",
        format!(
            "{}-{}",
            region.extended_start.get(),
            region.extended_end.get()
        ),
    );
    kv("java_padded", &java_padded);
    kv(
        "iterator",
        format!(
            "cigar={} start={} len={}",
            cigar_of(raw),
            raw.pos() + 1,
            bases_of(raw).len()
        ),
    );
    kv(
        "after_finalize",
        format!(
            "cigar={} start={} len={}",
            cigar_of(after_finalize),
            after_finalize.pos() + 1,
            bases_of(after_finalize).len()
        ),
    );
    kv(
        "rust_trim_padded",
        format!("{}-{}", trim.trim_start, trim.trim_end),
    );
    kv("java_variant_padded", &java_variant_padded);
    kv(
        "after_trim_hardclip",
        format!(
            "cigar={} start={} len={}",
            cigar_of(&clipped),
            clipped.pos() + 1,
            bases_of(&clipped).len()
        ),
    );

    assert_eq!(java_padded, "20:29455460-29455844");
    assert_eq!(
        (region.extended_start.get(), region.extended_end.get()),
        (29_455_460, 29_455_844)
    );
    assert_eq!(cigar_of(raw), "148M");
    assert_eq!(raw.pos() + 1, 29_455_509);
    assert_eq!(bases_of(raw).len(), 148);
    assert_eq!(bases_of(after_finalize), bases_of(raw));
    assert_eq!(cigar_of(after_finalize), "148M");
    assert_eq!(java_raw.cigar, "148M");
    assert_eq!(java_raw.start, 29_455_509);
    assert_eq!(java_raw.bases, bases_of(raw));
    assert_eq!(java_final.cigar, "148M");
    assert_eq!(java_final.start, 29_455_509);
    assert_eq!(java_final.bases, bases_of(raw));

    assert_eq!(java_variant_padded, "20:29455560-29455728");
    assert_eq!((trim.trim_start, trim.trim_end), (29_455_569, 29_455_724));
    assert_eq!(java_trim.cigar, "51H97M");
    assert_eq!(java_trim.start, 29_455_560);
    assert_eq!(java_trim.len, 97);
    assert_eq!(java_trim.end, 29_455_656);
    assert_eq!(cigar_of(&clipped), "60H88M");
    assert_eq!(clipped.pos() + 1, 29_455_569);
    assert_eq!(bases_of(&clipped).len(), 88);
    assert_eq!(&java_trim.bases[9..], bases_of(&clipped).as_slice());
    assert_eq!(&bases_of(raw)[51..60], EXTRA.as_bytes());
    assert_eq!(&java_trim.bases[..9], EXTRA.as_bytes());
    assert_eq!(&bases_of(raw)[51..], java_trim.bases.as_slice());
    assert_eq!(&bases_of(raw)[60..], bases_of(&clipped).as_slice());

    kv(
        "last_equal",
        "148M @ 29455509 after finalizeRegion/finalize_owned_bam_records",
    );
    kv(
        "first_divergent_op",
        "hardClipToRegion on the trimmed padded span",
    );
    kv(
        "java_op",
        "AssemblyRegion.trim:266 ReadClipper.hardClipToRegion(read, 29455560, 29455728) -> ReadClipper.hardClipByReferenceCoordinatesLeftTail(read, 29455559)",
    );
    kv(
        "rust_op",
        "clip_finalized_reads_in_place -> hard_clip_to_region(read, 29455569, 29455724) -> hard_clip_by_ref_left_tail(read, 29455568)",
    );
    kv("nine_bp", EXTRA);
    kv(
        "nine_bp_disposition",
        "Java retains raw[51..60]; Rust removes them",
    );
    kv("classification", "READ_HARDCLIP_TRIMMED_SPAN_DIVERGENCE");
    kv("production_change", "NONE");
    kv("counterfactual", "NOT_RUN");
}
