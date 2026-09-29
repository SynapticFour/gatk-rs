//! 6R.253: why Java trim span `20:29455560-29455728` vs Rust `20:29455569-29455724`
//! at ActiveFull `20:29455560-29455744` for `20:29455649 T/TGTTTG`.
//!
//! Actual GATK 4.4.0.0 executable dump (`hap-trim-at-loc` + `dumpTrimIntervalTrace`).
//! PRODUCTION CHANGE: NONE. PL ±1 is not investigated in this round.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r253_java_vs_rust_trim_span -- --nocapture --test-threads=1
//! HOLDOUT_6R253=1 cargo test -p gatk-haplotypecaller --test holdout_6r253_java_vs_rust_trim_span -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly_based_caller::assemble_reads_with_finalized;
use gatk_haplotypecaller::assembly_region_trimmer::{AssemblyRegionTrimmer, TrimVariant};
use gatk_haplotypecaller::reference_context::ReferenceContext;
use gatk_haplotypecaller::{
    begin_hap_list_observe, call_disposition, flatten_assembly_regions, take_hap_list_trim_span,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_455_649;
const A_INS: u64 = 29_455_644;
const JAVA_PADDED: (u64, u64) = (29_455_560, 29_455_728);
const RUST_PADDED: (u64, u64) = (29_455_569, 29_455_724);
const JAVA_TSV: &str = include_str!("forensic_6r253_java_trim_span.tsv");
const SNP_PAD: u64 = 20;
const INDEL_PAD: u64 = 75;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R253\t{key}\t{}", value.as_ref());
}

fn fields(line: &str) -> BTreeMap<String, String> {
    let mut out = BTreeMap::new();
    for p in line.split('\t').skip(2) {
        if let Some((k, v)) = p.split_once('=') {
            out.insert(k.to_string(), v.to_string());
        } else if p.contains(':') && !p.contains('=') {
            out.insert("span".to_string(), p.to_string());
        }
    }
    out
}

fn tsv_line(prefix: &str) -> &str {
    JAVA_TSV
        .lines()
        .find(|l| l.starts_with(prefix))
        .unwrap_or_else(|| panic!("missing {prefix}"))
}

fn pad_row(start: u64, alt: &str) -> BTreeMap<String, String> {
    for line in JAVA_TSV.lines() {
        if !line.contains("\tpad\t") {
            continue;
        }
        let f = fields(line);
        if f.get("start").map(String::as_str) == Some(&start.to_string())
            && f.get("alts").map(String::as_str) == Some(alt)
        {
            return f;
        }
    }
    panic!("missing Java pad row {start} {alt}");
}

fn replay_pad(events: &[TrimVariant], a_ins_pad: u64) -> (u64, u64) {
    let overlapping: Vec<&TrimVariant> = events
        .iter()
        .filter(|v| v.start <= 29_455_744 && v.end >= 29_455_560)
        .collect();
    let mut min_start = overlapping.iter().map(|v| v.start).min().unwrap();
    let mut max_end = overlapping.iter().map(|v| v.end).max().unwrap();
    for v in overlapping {
        let pad = if v.start == A_INS && v.is_indel {
            a_ins_pad
        } else if v.is_indel {
            INDEL_PAD
        } else {
            SNP_PAD
        };
        min_start = min_start.min(v.start.saturating_sub(pad).max(1));
        max_end = max_end.max(v.end.saturating_add(pad));
    }
    (min_start.max(29_455_460), max_end.min(29_455_844))
}

#[test]
fn forensic_6r253_java_oracle_trim_interval_objects() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv("pl_investigated", "NO");
    assert!(JAVA_TSV.contains("snpPadding=20"));
    assert!(JAVA_TSV.contains("indelPadding=75"));
    assert!(JAVA_TSV.contains("strPadding=75"));
    assert!(JAVA_TSV.contains("legacy=false"));
    assert_eq!(
        tsv_line("6R253\tactive\t"),
        "6R253\tactive\t20:29455560-29455744"
    );
    assert_eq!(
        tsv_line("6R253\tpadded\t"),
        "6R253\tpadded\t20:29455460-29455844"
    );
    assert_eq!(
        tsv_line("6R253\tref_window\t"),
        "6R253\tref_window\t20:29455460-29455844"
    );
    assert_eq!(
        tsv_line("6R253\tvariant_span_before_pad\t"),
        "6R253\tvariant_span_before_pad\t20:29455590-29455703"
    );
    assert_eq!(
        tsv_line("6R253\tpadded_variant_span\t"),
        "6R253\tpadded_variant_span\t20:29455560-29455728"
    );
    assert!(JAVA_TSV.contains("coordinate=reference"));
    assert_eq!(
        tsv_line("6R131\ttrim_span\t"),
        "6R131\ttrim_span\t20:29455560-29455728"
    );
    let a_at = pad_row(A_INS, "AT");
    assert_eq!(a_at["padding"], "84");
    assert!(a_at["str"].contains("unit=T"));
    assert!(a_at["str"].contains("longestSTR=9"));
    assert_eq!(a_at["start_minus_pad"], "29455560");
    assert_eq!(a_at["end_plus_pad"], "29455728");
    let tgtttg = pad_row(TARGET, "TGTTTG");
    assert_eq!(tgtttg["padding"], "75");
    assert_eq!(tgtttg["str"], "none");
    assert_eq!(tgtttg["end_plus_pad"], "29455724");
    kv(
        "java_first_span_setter",
        "A/AT 29455644 STR padding=84 (75+longestSTR=9)",
    );
}

#[test]
fn forensic_6r253_java_vs_rust_trim_span() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv("pl_investigated", "NO");
    kv(
        "classification_candidate",
        "TRIM_PADDING_SEMANTICS_DIVERGENCE",
    );

    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
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
        "rust_active",
        format!("20:{}-{}", region.start.get(), region.end.get()),
    );
    kv(
        "rust_padded_assembly",
        format!(
            "20:{}-{}",
            region.extended_start.get(),
            region.extended_end.get()
        ),
    );
    assert_eq!(region.start.get(), 29_455_560);
    assert_eq!(region.end.get(), 29_455_744);
    assert_eq!(region.extended_start.get(), 29_455_460);
    assert_eq!(region.extended_end.get(), 29_455_844);

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
    let trim_variants: Vec<TrimVariant> = assembled
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
    let overlapping: Vec<&TrimVariant> = trim_variants
        .iter()
        .filter(|v| v.overlaps_active_region(region))
        .collect();
    kv("rust_n_overlap", overlapping.len().to_string());
    let rust_has_a_ins = overlapping.iter().any(|v| v.start == A_INS && v.is_indel);
    let rust_has_tgtttg = overlapping.iter().any(|v| v.start == TARGET && v.is_indel);
    kv("rust_has_A_AT", rust_has_a_ins.to_string());
    kv("rust_has_TGTTTG", rust_has_tgtttg.to_string());
    assert!(rust_has_a_ins, "Rust EventMap union includes A/AT for trim");
    assert!(rust_has_tgtttg, "Rust EventMap union includes T/TGTTTG");
    let a_ins = overlapping
        .iter()
        .find(|v| v.start == A_INS && v.is_indel)
        .expect("A/AT");
    kv(
        "rust_A_AT_span",
        format!("{}-{} indel={}", a_ins.start, a_ins.end, a_ins.is_indel),
    );
    assert_eq!(a_ins.start, a_ins.end, "insertion VCF end == start");

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
    let rust_trim = trimmer.trim(region, &trim_variants, Some(&ref_ctx));
    let rust_ps = rust_trim.padded_variant_start.expect("padded start");
    let rust_pe = rust_trim.padded_variant_end.expect("padded end");
    kv("rust_trimmer_padded", format!("{rust_ps}-{rust_pe}"));
    kv(
        "rust_trimmer_variant",
        format!(
            "{}-{}",
            rust_trim.variant_start.unwrap(),
            rust_trim.variant_end.unwrap()
        ),
    );
    assert_eq!((rust_ps, rust_pe), RUST_PADDED);
    assert_eq!(rust_trim.variant_start, Some(29_455_590));
    assert_eq!(rust_trim.variant_end, Some(29_455_703));

    let rust_replay = replay_pad(&trim_variants, INDEL_PAD);
    let java_replay = replay_pad(&trim_variants, 84);
    kv(
        "rust_replay_indel75",
        format!("{}-{}", rust_replay.0, rust_replay.1),
    );
    kv(
        "java_replay_A_AT_pad84",
        format!("{}-{}", java_replay.0, java_replay.1),
    );
    assert_eq!(rust_replay, RUST_PADDED);
    assert_eq!(java_replay, JAVA_PADDED);

    begin_hap_list_observe();
    let _ = HaplotypeCallerEngine::call_region(region, &dict, &ref_fasta, &args)
        .expect("call")
        .expect("Some");
    let span = take_hap_list_trim_span().expect("trim span");
    kv(
        "production_call_region_trim",
        format!("{}-{}", span.trim_start, span.trim_end),
    );
    assert_eq!(span.trim_start, RUST_PADDED.0);
    assert_eq!(span.trim_end, RUST_PADDED.1);
    assert_ne!(span.trim_start, JAVA_PADDED.0);
    assert_ne!(span.trim_end, JAVA_PADDED.1);

    kv(
        "first_divergent_operation",
        "TandemRepeat.getNumTandemRepeatUnits on A/AT 29455644 → padding 84 vs Rust indel pad 75 (STR skipped: insertion start==end)",
    );
    kv("classification", "TRIM_PADDING_SEMANTICS_DIVERGENCE");
    kv("nine_left", "29455569-29455560=9 = Java longestSTR at A/AT");
    kv(
        "four_right",
        "29455728-29455724=4 = (A/AT+84) vs (TGTTTG+75)",
    );
    kv(
        "haplotype_length",
        "Java 169 ref + 5I = 174; Rust 156 ref + 5I = 161",
    );
}
