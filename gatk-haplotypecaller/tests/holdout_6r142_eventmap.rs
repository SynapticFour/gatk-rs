//! 6R.142: EventMap generation on the 6R.141-equivalent trimmed haplotype set.
//! Diagnostic-only `unbounded_diagnostic`. Production k-best policy unchanged.
//! Skipped unless `HOLDOUT_6R142=1`.
//!
//! ```text
//! HOLDOUT_6R142=1 GATK_RS_EXPERIMENTAL_KBEST_POLICY=unbounded_diagnostic \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r142_eventmap -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::event_map::{
    collect_variation_events, variation_events_for_haplotype, EventMap, VariationEvent,
};
use gatk_haplotypecaller::{
    assemble_reads_with_finalized, call_disposition, flatten_assembly_regions,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs, GenomePosition,
    Haplotype, ReadFilterParams, WalkerTraversalConfig,
};
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const JAVA_DUMP_REL: &str = "parity/giab/runs/local-pairhmm-diff/6r142_java.txt";
const TARGET: u64 = 29_456_196;
const JAVA_ACTIVE_START: u64 = 29_455_995;
const JAVA_ACTIVE_END: u64 = 29_456_293;
const JAVA_PAD_START: u64 = 29_455_895;
const JAVA_PAD_END: u64 = 29_456_393;
const JAVA_TRIM_START: u64 = 29_455_977;
const JAVA_TRIM_END: u64 = 29_456_301;
const EXPECTED_REF_HASH: &str = "659741f99b7f78a5";
const EXPECTED_TRIM_REF_HASH: &str = "9af99c7ab7cfc3ee";
const EXPECTED_FULL_REF_HASH: &str = "52a661e82a18ef87";
const EXPECTED_FINALIZED: usize = 379;
const JAVA_FULL_PAD_START: u64 = 29_455_395;
const MAX_MNP: usize = 0;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R142\t{key}\t{}", value.as_ref());
}

fn fnv1a64(bases: &[u8]) -> u64 {
    let mut h = 0xcbf29ce484222325u64;
    for &b in bases {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

fn hex(h: u64) -> String {
    format!("{h:016x}")
}

fn hap_hash(h: &Haplotype) -> String {
    hex(fnv1a64(&h.bases))
}

struct EnvGuard {
    key: &'static str,
    prior: Option<String>,
}

impl EnvGuard {
    fn set(key: &'static str, value: &str) -> Self {
        let prior = std::env::var(key).ok();
        std::env::set_var(key, value);
        Self { key, prior }
    }
}

impl Drop for EnvGuard {
    fn drop(&mut self) {
        match &self.prior {
            Some(v) => std::env::set_var(self.key, v),
            None => std::env::remove_var(self.key),
        }
    }
}

fn parse_kv_line(line: &str) -> BTreeMap<String, String> {
    let mut m = BTreeMap::new();
    for tok in line.split('\t') {
        if let Some((k, v)) = tok.split_once('=') {
            m.insert(k.to_string(), v.to_string());
        }
    }
    m
}

type EventKey = (u64, String, String);

fn event_key(start: u64, r: &str, a: &str) -> EventKey {
    (start, r.to_string(), a.to_string())
}

fn event_type(r: &str, a: &str) -> &'static str {
    if r.len() == 1 && a.len() == 1 {
        "SNP"
    } else if r.len() == a.len() {
        "MNP"
    } else {
        "INDEL"
    }
}

fn ve_key(e: &VariationEvent) -> EventKey {
    event_key(e.start_1based.get(), &e.ref_allele, &e.alt_allele)
}

struct JavaHap {
    #[allow(dead_code)]
    idx: usize,
    #[allow(dead_code)]
    hash: String,
    is_ref: bool,
    align: usize,
    cigar: String,
    events: Vec<EventKey>,
}

struct JavaDump {
    hap_order: Vec<String>,
    haps: BTreeMap<String, JavaHap>,
    union: BTreeSet<EventKey>,
    full_ref_hash: String,
    padded_ref_loc: String,
    max_mnp: usize,
    hap_n: usize,
}

fn load_java_trimmed(path: &Path) -> JavaDump {
    let text = std::fs::read_to_string(path).expect("6r142 java dump");
    let mut haps: BTreeMap<String, JavaHap> = BTreeMap::new();
    let mut hap_order = Vec::new();
    let mut union = BTreeSet::new();
    let mut full_ref_hash = String::new();
    let mut padded_ref_loc = String::new();
    let mut max_mnp = 0usize;
    let mut hap_n = 0usize;
    for line in text.lines() {
        if !line.starts_with("6R142\t") {
            continue;
        }
        let kind = line.split('\t').nth(1).unwrap_or("");
        let m = parse_kv_line(line);
        match kind {
            "input" if m.get("stage").map(|s| s.as_str()) == Some("trimmed") => {
                full_ref_hash = m.get("fullRef_hash").cloned().unwrap_or_default();
                padded_ref_loc = m.get("paddedRefLoc").cloned().unwrap_or_default();
                max_mnp = m
                    .get("maxMnpDistance")
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(0);
                hap_n = m.get("hap_n").and_then(|s| s.parse().ok()).unwrap_or(0);
            }
            "hap" if m.get("stage").map(|s| s.as_str()) == Some("trimmed") => {
                let hash = m.get("hash").cloned().expect("hash");
                hap_order.push(hash.clone());
                haps.insert(
                    hash.clone(),
                    JavaHap {
                        idx: m.get("idx").and_then(|s| s.parse().ok()).unwrap_or(0),
                        hash,
                        is_ref: m.get("isRef").map(|s| s == "true").unwrap_or(false),
                        align: m
                            .get("alignStart")
                            .and_then(|s| s.parse().ok())
                            .unwrap_or(0),
                        cigar: m.get("cigar").cloned().unwrap_or_else(|| ".".into()),
                        events: Vec::new(),
                    },
                );
            }
            "event" if m.get("stage").map(|s| s.as_str()) == Some("trimmed") => {
                let hash = m.get("hash").cloned().expect("hash");
                let start: u64 = m.get("start").and_then(|s| s.parse().ok()).unwrap_or(0);
                let r = m.get("ref").cloned().unwrap_or_default();
                let a = m.get("alt").cloned().unwrap_or_default();
                if let Some(h) = haps.get_mut(&hash) {
                    h.events.push(event_key(start, &r, &a));
                }
            }
            "union" if m.get("stage").map(|s| s.as_str()) == Some("trimmed") => {
                let start: u64 = m.get("start").and_then(|s| s.parse().ok()).unwrap_or(0);
                let r = m.get("ref").cloned().unwrap_or_default();
                let alts = m.get("alts").cloned().unwrap_or_default();
                for a in alts.split(',') {
                    if a != "." && !a.is_empty() {
                        union.insert(event_key(start, &r, a));
                    }
                }
            }
            _ => {}
        }
    }
    JavaDump {
        hap_order,
        haps,
        union,
        full_ref_hash,
        padded_ref_loc,
        max_mnp,
        hap_n,
    }
}

fn cigar_only_events(
    h: &Haplotype,
    ref_hap: &Haplotype,
    ref_bytes: &[u8],
    ref_loc: u64,
    contig: &str,
) -> Vec<VariationEvent> {
    EventMap::from_haplotype_and_reference(h, ref_hap, ref_bytes, ref_loc, MAX_MNP)
        .variation_events(contig, ref_loc)
}

fn set_of(events: &[VariationEvent]) -> BTreeSet<EventKey> {
    events.iter().map(ve_key).collect()
}

fn set_of_keys(keys: &[EventKey]) -> BTreeSet<EventKey> {
    keys.iter().cloned().collect()
}

fn dump_set_diff(label: &str, rust: &BTreeSet<EventKey>, java: &BTreeSet<EventKey>) {
    let common = rust.intersection(java).count();
    kv(
        "compare",
        format!(
            "label={label}\trust={}\tjava={}\tCOMMON={common}\tJAVA_ONLY={}\tRUST_ONLY={}",
            rust.len(),
            java.len(),
            java.difference(rust).count(),
            rust.difference(java).count()
        ),
    );
    for e in java.difference(rust).take(8) {
        kv(
            "java_only",
            format!(
                "label={label}\tstart={}\tref={}\talt={}\ttype={}",
                e.0,
                e.1,
                e.2,
                event_type(&e.1, &e.2)
            ),
        );
    }
    for e in rust.difference(java).take(8) {
        kv(
            "rust_only",
            format!(
                "label={label}\tstart={}\tref={}\talt={}\ttype={}",
                e.0,
                e.1,
                e.2,
                event_type(&e.1, &e.2)
            ),
        );
    }
}

#[test]
fn holdout_6r142_eventmap() {
    if std::env::var("HOLDOUT_6R142").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R142=1");
        return;
    }
    let _policy = EnvGuard::set("GATK_RS_EXPERIMENTAL_KBEST_POLICY", "unbounded_diagnostic");
    kv(
        "parity_mode",
        "unbounded_diagnostic\tnote=diagnostic_not_production_default",
    );

    let root = repo_root();
    let java = load_java_trimmed(&root.join(JAVA_DUMP_REL));
    kv("java_trimmed_haps", java.hap_n.to_string());
    kv("java_union_n", java.union.len().to_string());
    kv("java_fullRef_hash", java.full_ref_hash.as_str());
    kv("java_paddedRefLoc", java.padded_ref_loc.as_str());
    kv("java_maxMnp", java.max_mnp.to_string());
    assert_eq!(java.hap_n, 25);
    assert_eq!(java.hap_order.len(), 25);
    assert_eq!(java.full_ref_hash, EXPECTED_FULL_REF_HASH);
    assert_eq!(java.max_mnp, 0);

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
    let covering = regions
        .iter()
        .find(|r| {
            matches!(
                call_disposition(r),
                AssemblyRegionCallDisposition::ActiveFull
            ) && r.start.get() <= TARGET
                && r.end.get() >= TARGET
        })
        .expect("ActiveFull")
        .clone();
    let mut java_bounds = covering;
    java_bounds.start = GenomePosition::new_1based(JAVA_ACTIVE_START);
    java_bounds.end = GenomePosition::new_1based(JAVA_ACTIVE_END);
    java_bounds.extended_start = GenomePosition::new_1based(JAVA_PAD_START);
    java_bounds.extended_end = GenomePosition::new_1based(JAVA_PAD_END);

    let args = CallRegionArgs::strict_java();
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let mut owned = java_bounds.clone();
    let assembled =
        assemble_reads_with_finalized(&mut owned, &dict, &mut ref_cache, &args.assemble)
            .expect("assemble");
    assert_eq!(
        hex(fnv1a64(assembled.assembly.reference_bases())),
        EXPECTED_REF_HASH
    );
    assert_eq!(assembled.finalized_reads.len(), EXPECTED_FINALIZED);

    let full_ref_read = gatk_haplotypecaller::assembly_region_finalize::assembly_reference_read(
        &dict,
        &mut ref_cache,
        &java_bounds,
    )
    .expect("full pad ref");
    kv(
        "rust_fullRef",
        format!(
            "len={}\thash={}\tpad_start_expected={}",
            full_ref_read.bases.len(),
            hex(fnv1a64(&full_ref_read.bases)),
            JAVA_FULL_PAD_START
        ),
    );
    kv(
        "rust_graphRef",
        format!(
            "len={}\thash={}\tpad_start={}",
            assembled.assembly.reference_bases().len(),
            hex(fnv1a64(assembled.assembly.reference_bases())),
            assembled.assembly.padded_reference_start_1based()
        ),
    );

    let mut trim_region = java_bounds.clone();
    trim_region.extended_start = GenomePosition::new_1based(JAVA_TRIM_START);
    trim_region.extended_end = GenomePosition::new_1based(JAVA_TRIM_END);
    let trimmed = assembled.assembly.trim_to(&trim_region).expect("trim_to");
    let rust_order: Vec<String> = trimmed.haplotypes.iter().map(hap_hash).collect();
    kv(
        "hap_order",
        format!(
            "rust_n={}\tjava_n={}\tidentical={}",
            rust_order.len(),
            java.hap_order.len(),
            rust_order == java.hap_order
        ),
    );
    assert_eq!(rust_order.len(), 25);
    assert_eq!(rust_order, java.hap_order, "6R.141 trimmed order invariant");

    let contig = "20";
    let (graph_ref, graph_pad) = trimmed.event_map_reference();
    kv(
        "eventmap_ref_input",
        format!(
            "graph_len={}\tgraph_hash={}\tgraph_pad={}\tmax_mnp={}\tstrict_java={}",
            graph_ref.len(),
            hex(fnv1a64(graph_ref)),
            graph_pad,
            trimmed.max_mnp_distance(),
            true
        ),
    );

    let ref_hap = trimmed
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .expect("ref hap");
    kv(
        "rust_refHap",
        format!(
            "hash={}\tlen={}\talign={}\tcigar={}",
            hap_hash(ref_hap),
            ref_hap.bases.len(),
            ref_hap.alignment_start_hap_wrt_ref,
            ref_hap
                .cigar
                .as_ref()
                .map(|c| c.to_gatk_string())
                .unwrap_or_else(|| ".".into())
        ),
    );
    assert_eq!(hap_hash(ref_hap), EXPECTED_TRIM_REF_HASH);

    let mut cigar_vs_java_mismatch_haps = 0usize;
    let mut prod_vs_java_mismatch_haps = 0usize;
    let mut prod_vs_cigar_mismatch_haps = 0usize;
    let mut align_mismatch = 0usize;
    let mut cigar_mismatch = 0usize;
    let mut isref_mismatch = 0usize;
    let mut first_div_hap = String::from("none");
    let mut first_div_event = String::from("none");
    let mut cigar_union = BTreeSet::new();
    let mut prod_union = BTreeSet::new();

    for (i, h) in trimmed.haplotypes.iter().enumerate() {
        let hash = hap_hash(h);
        let jh = java.haps.get(&hash).expect("java hap");
        let rust_cigar = h
            .cigar
            .as_ref()
            .map(|c| c.to_gatk_string())
            .unwrap_or_else(|| ".".into());
        if rust_cigar != jh.cigar {
            cigar_mismatch += 1;
        }
        if h.is_reference != jh.is_ref {
            isref_mismatch += 1;
        }
        if h.alignment_start_hap_wrt_ref != jh.align {
            align_mismatch += 1;
        }

        let cigar_evs = cigar_only_events(h, ref_hap, graph_ref, graph_pad, contig);
        let prod_evs =
            variation_events_for_haplotype(h, ref_hap, graph_ref, graph_pad, MAX_MNP, contig);
        let cigar_set = set_of(&cigar_evs);
        let prod_set = set_of(&prod_evs);
        let java_set = set_of_keys(&jh.events);
        cigar_union.extend(cigar_set.iter().cloned());
        prod_union.extend(prod_set.iter().cloned());

        kv(
            "hap",
            format!(
                "idx={i}\thash={hash}\tisRef={}\talign={}\tjava_align={}\tcigar={rust_cigar}\tcigar_n={}\tprod_n={}\tjava_n={}",
                h.is_reference,
                h.alignment_start_hap_wrt_ref,
                jh.align,
                cigar_evs.len(),
                prod_evs.len(),
                jh.events.len()
            ),
        );
        for (evi, e) in cigar_evs.iter().enumerate() {
            kv(
                "event",
                format!(
                    "kind=cigar_only\thash={hash}\tevi={evi}\tstart={}\tend={}\tref={}\talt={}\ttype={}",
                    e.start_1based.get(),
                    e.end_1based.get(),
                    e.ref_allele,
                    e.alt_allele,
                    event_type(&e.ref_allele, &e.alt_allele)
                ),
            );
        }

        if cigar_set != java_set {
            cigar_vs_java_mismatch_haps += 1;
            if first_div_hap == "none" {
                first_div_hap = hash.clone();
                if let Some(e) = java_set.difference(&cigar_set).next() {
                    first_div_event = format!("JAVA_ONLY {}:{}>{}", e.0, e.1, e.2);
                } else if let Some(e) = cigar_set.difference(&java_set).next() {
                    first_div_event = format!("RUST_ONLY {}:{}>{}", e.0, e.1, e.2);
                }
            }
            dump_set_diff(&format!("hap_{hash}_cigar"), &cigar_set, &java_set);
        }
        if prod_set != java_set {
            prod_vs_java_mismatch_haps += 1;
            dump_set_diff(&format!("hap_{hash}_prod"), &prod_set, &java_set);
        }
        if prod_set != cigar_set {
            prod_vs_cigar_mismatch_haps += 1;
            dump_set_diff(&format!("hap_{hash}_prod_vs_cigar"), &prod_set, &cigar_set);
        }
    }

    let stored_union = set_of(trimmed.variation_events());
    let collect_union = set_of(&collect_variation_events(
        &trimmed.haplotypes,
        graph_ref,
        graph_pad,
        contig,
        MAX_MNP,
    ));
    dump_set_diff("union_cigar_vs_java", &cigar_union, &java.union);
    dump_set_diff("union_prod_vs_java", &prod_union, &java.union);
    dump_set_diff("union_stored_vs_java", &stored_union, &java.union);
    dump_set_diff("union_collect_vs_java", &collect_union, &java.union);

    // Counterfactual: Java fullReferenceWithPadding + Java-style alignStart (graph offset + 500).
    let mut fullpad_union = BTreeSet::new();
    let mut fullpad_mismatch = 0usize;
    let owned_full_ref;
    let full_ref_hap = {
        owned_full_ref = {
            let mut rh = Haplotype::new(full_ref_read.bases.clone(), true);
            rh.alignment_start_hap_wrt_ref = 500;
            rh
        };
        &owned_full_ref
    };
    for h in &trimmed.haplotypes {
        let mut lifted = h.clone();
        lifted.alignment_start_hap_wrt_ref = h.alignment_start_hap_wrt_ref.saturating_add(500);
        let evs = cigar_only_events(
            &lifted,
            full_ref_hap,
            &full_ref_read.bases,
            JAVA_FULL_PAD_START,
            contig,
        );
        let set = set_of(&evs);
        fullpad_union.extend(set.iter().cloned());
        let java_set = set_of_keys(&java.haps.get(&hap_hash(h)).expect("jh").events);
        if set != java_set {
            fullpad_mismatch += 1;
        }
    }
    dump_set_diff("union_fullpad_lifted_vs_java", &fullpad_union, &java.union);

    let hap_maps_equal = cigar_vs_java_mismatch_haps == 0 && prod_vs_java_mismatch_haps == 0;
    let union_equal = cigar_union == java.union && prod_union == java.union;
    let (classification, first_divergence) = if rust_order != java.hap_order {
        ("EVENTMAP_GENERATION_DIVERGENCE", "trimmed haplotype order")
    } else if cigar_vs_java_mismatch_haps > 0 {
        (
            "EVENTMAP_GENERATION_DIVERGENCE",
            "processCigarForInitialEvents analogue (EventMap::from_haplotype_and_reference)",
        )
    } else if prod_vs_cigar_mismatch_haps > 0 {
        (
            "EVENTMAP_GENERATION_DIVERGENCE",
            "variation_events_for_haplotype extras vs CIGAR EventMap",
        )
    } else if !union_equal {
        (
            "EVENTMAP_GENERATION_DIVERGENCE",
            "EventMap union / getVariationEvents",
        )
    } else if align_mismatch > 0 {
        (
            "EVENTMAP_REPRESENTATION_ONLY",
            "alignmentStartHapwrtRef (500-offset vs graph-pad 0-offset; genomic events identical)",
        )
    } else {
        ("NO_DIVERGENCE", "none")
    };

    kv(
        "classification",
        format!(
            "class={classification}\tfirst_divergence={first_divergence}\thap_maps_equal={hap_maps_equal}\tunion_equal={union_equal}\tcigar_vs_java_haps={cigar_vs_java_mismatch_haps}\tprod_vs_java_haps={prod_vs_java_mismatch_haps}\tprod_vs_cigar_haps={prod_vs_cigar_mismatch_haps}\talign_mismatch={align_mismatch}\tcigar_mismatch={cigar_mismatch}\tisref_mismatch={isref_mismatch}\tfullpad_mismatch_haps={fullpad_mismatch}\tfirst_div_hap={first_div_hap}\tfirst_div_event={first_div_event}"
        ),
    );
    kv(
        "summary",
        format!(
            "haplotypes_java={}/25\thaplotypes_rust={}/25\teventmaps_equal={}\tunion_java={}\tunion_rust_cigar={}",
            java.hap_n,
            rust_order.len(),
            hap_maps_equal && union_equal,
            java.union.len(),
            cigar_union.len()
        ),
    );

    assert_eq!(cigar_mismatch, 0);
    assert_eq!(isref_mismatch, 0);
    assert_eq!(cigar_vs_java_mismatch_haps, 0);
    assert_eq!(prod_vs_java_mismatch_haps, 0);
    assert_eq!(prod_vs_cigar_mismatch_haps, 0);
    assert_eq!(cigar_union, java.union);
    assert_eq!(prod_union, java.union);
    assert_eq!(collect_union, java.union);
}
