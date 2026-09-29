//! 6R.232: first operation that creates/retains Java variation event
//! `20:29455314` while Rust does not. Target site remains `20:29455379 G/A`.
//!
//! Proof-only. Frozen Java 4.4.0.0 SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! PRODUCTION CHANGE: NONE.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r232_missing_variation_event_boundary -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::assembly_based_caller::assemble_reads_with_finalized;
use gatk_haplotypecaller::assembly_region_trimmer::{AssemblyRegionTrimmer, TrimVariant};
use gatk_haplotypecaller::cigar::CigarOperator;
use gatk_haplotypecaller::event_map::{
    build_per_haplotype_variation_events, collect_variation_events,
};
use gatk_haplotypecaller::haplotype::Haplotype;
use gatk_haplotypecaller::read_unclip::alignment_end_1based;
use gatk_haplotypecaller::reference_context::ReferenceContext;
use gatk_haplotypecaller::{
    begin_hap_list_observe, call_disposition, flatten_assembly_regions, take_hap_list_trim_span,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use rust_htslib::bam::record::Cigar;
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

const JAVA_PIN: &str = "4.4.0.0 / 2dbc025821bc5f686c423ff332a41e6cef892a77";
const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_455_379;
const EVENT_POS: u64 = 29_455_314;
const EVENT_REF: &str = "G";
const EVENT_ALT: &str = "C";
const JAVA_CARRIER_UNTRIM: &str = "c7acc50dfb9f9ecc";
const SNP_PADDING: u64 = 20;
const JAVA_PADDED_START: u64 = 29_455_294;
const RUST_PADDED_START: u64 = 29_455_355;
const JAVA_TSV: &str = include_str!("forensic_6r232_java_eventmap.tsv");

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    eprintln!("6R232\t{key}\t{}", value.as_ref());
}

fn fnv1a64_hex(data: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in data {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{h:016x}")
}

fn parse_kv_fields(line: &str) -> BTreeMap<String, String> {
    let mut out = BTreeMap::new();
    for p in line.split('\t').skip(2) {
        if let Some((k, v)) = p.split_once('=') {
            out.insert(k.to_string(), v.to_string());
        }
    }
    out
}

fn hap_base_at(h: &Haplotype, pos: u64) -> Option<u8> {
    let gl = h.genome_loc?;
    if pos < gl.start_1based() || pos > gl.end_1based() {
        return None;
    }
    let cigar = h.cigar.as_ref()?;
    let mut ref_pos = gl.start_1based();
    let mut hap_pos = 0usize;
    for el in &cigar.elements {
        match el.operator {
            CigarOperator::Match => {
                for _ in 0..el.length {
                    if ref_pos == pos {
                        return h.bases.get(hap_pos).copied();
                    }
                    ref_pos += 1;
                    hap_pos += 1;
                }
            }
            CigarOperator::Insertion | CigarOperator::SoftClip => {
                hap_pos += el.length;
            }
            CigarOperator::Deletion => {
                for _ in 0..el.length {
                    if ref_pos == pos {
                        return None;
                    }
                    ref_pos += 1;
                }
            }
            CigarOperator::HardClip => {}
        }
    }
    None
}

fn pileup_at(
    reads: &[gatk_haplotypecaller::shared_bam::SharedBamRecord],
    pos: u64,
) -> BTreeMap<u8, usize> {
    let mut counts: BTreeMap<u8, usize> = BTreeMap::new();
    for rec in reads {
        let r = rec.as_ref();
        if r.tid() < 0 || r.is_unmapped() {
            continue;
        }
        let mut ref_pos = (r.pos() + 1) as u64;
        let seq = r.seq();
        let mut read_pos = 0usize;
        for c in r.cigar().iter() {
            match c {
                Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                    for _ in 0..*n {
                        if ref_pos == pos {
                            *counts.entry(seq.as_bytes()[read_pos]).or_insert(0) += 1;
                        }
                        ref_pos += 1;
                        read_pos += 1;
                    }
                }
                Cigar::Ins(n) | Cigar::SoftClip(n) => read_pos += *n as usize,
                Cigar::Del(n) | Cigar::RefSkip(n) => ref_pos += u64::from(*n),
                Cigar::HardClip(_) | Cigar::Pad(_) => {}
            }
            if ref_pos > pos {
                break;
            }
        }
        let _ = alignment_end_1based(r);
    }
    counts
}

#[test]
fn forensic_6r232_missing_variation_event_boundary() {
    kv("java_pin", JAVA_PIN);
    kv("production_change", "NONE");
    kv("target", "20:29455379 G/A");
    kv(
        "predecessor",
        "6R.231 READ_INTERVAL_INPUT_DIVERGENCE CLOSED",
    );
    kv(
        "java_source",
        "EventMap.buildEventMapsForHaplotypes / AssemblyResultSet.getVariationEvents (SHA 2dbc0258)",
    );

    let mut java_union: Vec<(u64, u64, String, String)> = Vec::new();
    let mut saw_carrier = false;
    for line in JAVA_TSV.lines() {
        let f = parse_kv_fields(line);
        if line.contains("\tinput\t") && f.get("stage").map(String::as_str) == Some("untrimmed") {
            kv("java_untrimmed_hap_n", f.get("hap_n").unwrap().as_str());
            kv(
                "java_untrimmed_region",
                format!(
                    "{} padded={}",
                    f.get("region").unwrap(),
                    f.get("region_padded").unwrap()
                ),
            );
        }
        if line.contains("\tunion\t") && f.get("stage").map(String::as_str) == Some("untrimmed") {
            java_union.push((
                f["start"].parse().unwrap(),
                f["end"].parse().unwrap(),
                f["ref"].clone(),
                f["alts"].clone(),
            ));
        }
        if line.contains("\tevent\t")
            && f.get("start").map(String::as_str) == Some("29455314")
            && f.get("stage").map(String::as_str) == Some("untrimmed")
        {
            saw_carrier = true;
            assert_eq!(f["ref"], EVENT_REF);
            assert_eq!(f["alt"], EVENT_ALT);
            assert_eq!(f["type"], "SNP");
            assert_eq!(f["hash"], JAVA_CARRIER_UNTRIM);
            kv(
                "java_event_identity",
                format!(
                    "20:{EVENT_POS} {EVENT_REF}>{EVENT_ALT} SNP hap_hash={} isRef=false cigar=460M",
                    f["hash"]
                ),
            );
        }
        if line.contains("\tvariant_span\t") {
            kv("java_variant_span", line.split('\t').last().unwrap());
        }
        if line.contains("\ttrim_span\t") {
            kv("java_trim_span", line.split('\t').last().unwrap());
        }
    }
    assert!(
        saw_carrier,
        "Java fixture must pin 29455314 G>C on carrier hap"
    );
    let java_at_event: Vec<_> = java_union
        .iter()
        .filter(|(s, _, _, _)| *s == EVENT_POS)
        .cloned()
        .collect();
    assert_eq!(java_at_event.len(), 1);
    assert_eq!(java_at_event[0].2, EVENT_REF);
    assert_eq!(java_at_event[0].3, EVENT_ALT);
    kv(
        "java_event_source",
        "single alternate haplotype (untrimmed idx=115, 460M vs ref); reference haplotype event_n=0; CIGAR-derived EventMap SNPs from M-mismatch bases, not later getAllVariantContexts synthesis",
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

    let pile = pileup_at(&region.reads, EVENT_POS);
    kv(
        "read_pileup_29455314",
        pile.iter()
            .map(|(b, n)| format!("{}:{n}", *b as char))
            .collect::<Vec<_>>()
            .join(","),
    );
    let mut c_qnames = Vec::new();
    for rec in &region.reads {
        let r = rec.as_ref();
        if r.tid() < 0 || r.is_unmapped() {
            continue;
        }
        let mut ref_pos = (r.pos() + 1) as u64;
        let seq = r.seq();
        let mut read_pos = 0usize;
        let mut saw_c = false;
        for c in r.cigar().iter() {
            match c {
                Cigar::Match(n) | Cigar::Equal(n) | Cigar::Diff(n) => {
                    for _ in 0..*n {
                        if ref_pos == EVENT_POS && seq.as_bytes()[read_pos] == b'C' {
                            saw_c = true;
                        }
                        ref_pos += 1;
                        read_pos += 1;
                    }
                }
                Cigar::Ins(n) | Cigar::SoftClip(n) => read_pos += *n as usize,
                Cigar::Del(n) | Cigar::RefSkip(n) => ref_pos += u64::from(*n),
                Cigar::HardClip(_) | Cigar::Pad(_) => {}
            }
        }
        if saw_c {
            c_qnames.push(format!(
                "{}:{}",
                String::from_utf8_lossy(r.qname()),
                r.flags()
            ));
        }
    }
    kv("read_C_support_qnames", c_qnames.join(","));
    assert_eq!(c_qnames.len(), 4);

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
    let haps = &assembled.assembly.haplotypes;
    let (full_ref, full_pad) = assembled.assembly.event_map_reference();
    let max_mnp = assembled.assembly.max_mnp_distance();
    kv("rust_untrimmed_hap_n", haps.len().to_string());
    kv("rust_max_mnp_distance", max_mnp.to_string());
    kv("rust_eventmap_pad_start", full_pad.to_string());

    let rust_hashes: BTreeSet<String> = haps.iter().map(|h| fnv1a64_hex(&h.bases)).collect();
    kv(
        "rust_has_java_carrier_hash",
        rust_hashes.contains(JAVA_CARRIER_UNTRIM).to_string(),
    );
    assert!(
        rust_hashes.contains(JAVA_CARRIER_UNTRIM),
        "6R.239: K=128 now includes Java carrier {JAVA_CARRIER_UNTRIM}"
    );

    let mut hap_base_counts: BTreeMap<String, usize> = BTreeMap::new();
    let mut n_cover = 0usize;
    let mut n_c = 0usize;
    for h in haps {
        if let Some(b) = hap_base_at(h, EVENT_POS) {
            n_cover += 1;
            *hap_base_counts.entry((b as char).to_string()).or_insert(0) += 1;
            if b == b'C' {
                n_c += 1;
                kv(
                    "rust_hap_with_C",
                    format!(
                        "hash={} isRef={} len={}",
                        fnv1a64_hex(&h.bases),
                        h.is_reference,
                        h.bases.len()
                    ),
                );
            }
        }
    }
    kv("rust_haps_covering_29455314", n_cover.to_string());
    kv(
        "rust_hap_bases_at_29455314",
        hap_base_counts
            .iter()
            .map(|(b, n)| format!("{b}:{n}"))
            .collect::<Vec<_>>()
            .join(","),
    );
    kv("rust_haps_with_alt_C", n_c.to_string());
    assert!(n_c >= 1, "6R.239: carrier in K=128 carries C at 29455314");

    let per_hap = build_per_haplotype_variation_events(
        haps,
        full_ref,
        full_pad,
        max_mnp,
        region.contig.as_str(),
    );
    let mut per_hap_n = 0usize;
    for i in 0..per_hap.hap_count() {
        for e in per_hap.events_for(i) {
            if e.start_1based.get() == EVENT_POS {
                per_hap_n += 1;
                kv(
                    "rust_per_hap_event_314",
                    format!(
                        "hap={i} {}>{} isRef={}",
                        e.ref_allele, e.alt_allele, haps[i].is_reference
                    ),
                );
            }
        }
    }
    kv("rust_per_hap_eventmap_314_n", per_hap_n.to_string());
    assert!(
        per_hap_n >= 1,
        "6R.239: EventMap emits 29455314 from the carrier"
    );

    let union = collect_variation_events(haps, full_ref, full_pad, region.contig.as_str(), max_mnp);
    let stored = assembled.assembly.variation_events();
    let union_has = union.iter().any(|e| e.start_1based.get() == EVENT_POS);
    let stored_has = stored.iter().any(|e| e.start_1based.get() == EVENT_POS);
    kv("rust_union_contains_314", union_has.to_string());
    kv(
        "rust_stored_variation_events_contains_314",
        stored_has.to_string(),
    );
    assert!(union_has, "6R.239: union includes 29455314 G>C");
    assert!(
        stored_has,
        "6R.239: stored variation events include 29455314 G>C"
    );

    kv(
        "rust_spine_note",
        "production call_region includes discover_parity_spine_snp_events in trim_variants; matching padded start 29455355 proves spine did not insert 29455314",
    );

    for pos in [
        29_455_296, EVENT_POS, 29_455_328, 29_455_337, 29_455_375, TARGET,
    ] {
        let java = java_union
            .iter()
            .find(|(s, _, _, _)| *s == pos)
            .map(|(_, _, r, a)| format!("{r}>{a}"));
        let rust = stored
            .iter()
            .filter(|e| e.start_1based.get() == pos)
            .map(|e| format!("{}>{}", e.ref_allele, e.alt_allele))
            .collect::<Vec<_>>();
        let rust_hap_n = (0..per_hap.hap_count())
            .filter(|i| {
                per_hap
                    .events_for(*i)
                    .iter()
                    .any(|e| e.start_1based.get() == pos)
            })
            .count();
        let java_hap_n =
            if pos == EVENT_POS || pos == 29_455_296 || pos == 29_455_328 || pos == 29_455_337 {
                1
            } else {
                usize::MAX
            };
        kv(
            "window_event",
            format!(
                "pos={pos}\tjava={}\trust={}\trust_hap_n={rust_hap_n}\tjava_hap_n={}",
                java.as_deref().unwrap_or("absent"),
                if rust.is_empty() {
                    "absent".to_string()
                } else {
                    rust.join(",")
                },
                if java_hap_n == usize::MAX {
                    "union".to_string()
                } else {
                    java_hap_n.to_string()
                }
            ),
        );
    }

    let trim_variants: Vec<TrimVariant> = stored
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
    let rust_leftmost = overlapping.iter().map(|v| v.start).min().unwrap();
    kv("rust_leftmost_trim_event", rust_leftmost.to_string());
    assert_eq!(rust_leftmost, EVENT_POS);
    assert!(
        overlapping.iter().any(|v| v.start == EVENT_POS),
        "6R.239: trim-span sees 29455314 G>C"
    );

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
    kv(
        "rust_padded",
        format!(
            "{}-{}",
            rust_trim.padded_variant_start.unwrap(),
            rust_trim.padded_variant_end.unwrap()
        ),
    );
    assert_eq!(rust_trim.padded_variant_start.unwrap(), JAVA_PADDED_START);
    assert_eq!(
        rust_trim.padded_variant_start.unwrap(),
        rust_leftmost - SNP_PADDING
    );

    begin_hap_list_observe();
    let _ = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let span = take_hap_list_trim_span().expect("trim span");
    kv(
        "production_trim_span",
        format!("{}-{}", span.trim_start, span.trim_end),
    );
    assert_eq!(span.trim_start, JAVA_PADDED_START);

    let mut rust_plus_314 = trim_variants.clone();
    rust_plus_314.push(TrimVariant {
        contig: region.contig.clone(),
        start: EVENT_POS,
        end: EVENT_POS,
        is_indel: false,
        ref_allele: String::new(),
        alt_allele: String::new(),
    });
    let plus = trimmer.trim(region, &rust_plus_314, Some(&ref_ctx));
    kv(
        "diagnostic_rust_plus_314_padded_start",
        plus.padded_variant_start.unwrap().to_string(),
    );
    assert_eq!(plus.padded_variant_start.unwrap(), JAVA_PADDED_START);

    let java_overlapping: Vec<TrimVariant> = java_union
        .iter()
        .filter(|(s, e, _, _)| *s <= region.end.get() && *e >= region.start.get())
        .map(|(s, e, r, a)| TrimVariant {
            contig: region.contig.clone(),
            start: *s,
            end: *e,
            is_indel: r.len() != a.len() || s != e,
            ref_allele: r.clone(),
            alt_allele: a.clone(),
        })
        .collect();
    let java_full = trimmer.trim(region, &java_overlapping, Some(&ref_ctx));
    kv(
        "diagnostic_java_event_set_padded_start",
        java_full.padded_variant_start.unwrap().to_string(),
    );
    assert_eq!(java_full.padded_variant_start.unwrap(), JAVA_PADDED_START);

    let java_minus_314: Vec<TrimVariant> = java_overlapping
        .iter()
        .filter(|v| v.start != EVENT_POS)
        .cloned()
        .collect();
    let minus = trimmer.trim(region, &java_minus_314, Some(&ref_ctx));
    kv(
        "diagnostic_java_minus_314_padded_start",
        minus.padded_variant_start.unwrap().to_string(),
    );
    assert_eq!(minus.padded_variant_start.unwrap(), 29_455_308);
    assert_ne!(
        minus.padded_variant_start.unwrap(),
        RUST_PADDED_START,
        "removing only 29455314 does not make Java match Rust; 29455328/29455337 remain"
    );

    let java_minus_left_cluster: Vec<TrimVariant> = java_overlapping
        .iter()
        .filter(|v| v.start >= 29_455_375)
        .cloned()
        .collect();
    let minus_cluster = trimmer.trim(region, &java_minus_left_cluster, Some(&ref_ctx));
    kv(
        "diagnostic_java_minus_314_328_337_padded_start",
        minus_cluster.padded_variant_start.unwrap().to_string(),
    );
    assert_eq!(
        minus_cluster.padded_variant_start.unwrap(),
        RUST_PADDED_START
    );

    kv(
        "case",
        "6R.239 CLOSED — Java 29455314 G>C is now in Rust K=128 carrier EventMap; trim span matches Java 29455294",
    );
    kv(
        "stage_haplotype_bases",
        "Java alt hap C at 29455314; Rust all covering haps lack C (first divergence)",
    );
    kv("stage_per_hap_eventmap", "Java 1/128; Rust 0");
    kv("stage_eventmap_union", "Java present G>C; Rust absent");
    kv(
        "stage_merged_vc",
        "same as union for trim (getVariationEvents)",
    );
    kv(
        "stage_variant_span",
        "Java 29455314-29455523; Rust 29455375-29455523",
    );
    kv(
        "stage_padded_interval",
        "Java 29455294-29455584; Rust 29455355-29455581",
    );
    kv(
        "interval_arithmetic",
        "SNP pad 20 identical; left split is event membership, not formula",
    );
    kv("classification", "HAPLOTYPE_CONSTRUCTION_DIVERGENCE");
    kv(
        "first_causal_arrow",
        "Java assembly materializes haplotype c7acc50dfb9f9ecc (460M) whose EventMap emits 20:29455314 G>C; Rust never materializes a haplotype with C at that coordinate, so EventMap/union/trim never see the event. Do not patch EventMap.",
    );
    kv(
        "clipping_downstream",
        "6R.231 clip interval is the padded leftmost event; padding formula is shared",
    );
    kv(
        "pairhmm_downstream",
        "6R.230 read-sequence split is the clip of that interval",
    );
    kv(
        "haplotype_split_60_vs_54",
        "concurrent observation on the trimmed list; this round only proves the Java-only carrier of 29455314 is absent from Rust. EventMap is not the first miss.",
    );
    kv("production_src", "NONE");
}
