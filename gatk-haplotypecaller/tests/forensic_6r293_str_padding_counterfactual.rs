//! 6R.293: coordinate-only counterfactual. Add Java's proven STR contribution
//! of 9 to the single event `A>AT` at 20:29455644 and replay the trim loop.
//! The replay still uses the pre-6R.310 helper. 6R.310 made the live
//! production span the span this counterfactual reconstructed.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r293_str_padding_counterfactual -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary};
use gatk_haplotypecaller::{
    begin_trim_calc_observe, call_disposition, flatten_assembly_regions,
    longest_str_len_at_variant, take_trim_calc_snap, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams,
    ReferenceContext, TrimCalcSnap, WalkerTraversalConfig,
};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_455_649;
const SNP_PAD: u64 = 20;
const INDEL_PAD: u64 = 75;
const JAVA_STR: u64 = 9;
const JAVA_SPAN: (u64, u64) = (29_455_560, 29_455_728);
/// Span the pre-6R.310 helper produced. The replay below still reconstructs it.
const PRE_6R310_SPAN: (u64, u64) = (29_455_569, 29_455_724);

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R293\t{key}\t{}", value.as_ref());
}

fn span(bounds: (u64, u64)) -> String {
    format!("{}-{}", bounds.0, bounds.1)
}

#[derive(Clone)]
struct Event {
    start: u64,
    end: u64,
    is_indel: bool,
    ref_al: String,
    alt_al: String,
}

fn is_target(ev: &Event) -> bool {
    ev.start == 29_455_644
        && ev.end == 29_455_644
        && ev.ref_al == "A"
        && ev.alt_al == "AT"
        && ev.is_indel
}

fn is_next_indel(ev: &Event) -> bool {
    ev.start == 29_455_649
        && ev.end == 29_455_649
        && ev.ref_al == "T"
        && ev.alt_al == "TGTTTG"
        && ev.is_indel
}

/// Overlapping trim inputs in the order `trim_modern` walks them.
fn overlapping_in_trim_order(snap: &TrimCalcSnap) -> Vec<Event> {
    let mut events = snap.events.clone();
    let mut out = Vec::new();
    for v in &snap.trim_vars {
        if v.start > snap.active_end || v.end < snap.active_start {
            continue;
        }
        let idx = events
            .iter()
            .position(|e| e.start == v.start && e.end == v.end && e.is_indel == v.is_indel);
        let (ref_al, alt_al) = if let Some(i) = idx {
            let e = events.remove(i);
            (e.ref_al, e.alt_al)
        } else {
            (String::new(), String::new())
        };
        out.push(Event {
            start: v.start,
            end: v.end,
            is_indel: v.is_indel,
            ref_al,
            alt_al,
        });
    }
    out
}

fn padding(ev: &Event, ctx: &ReferenceContext, inject_target: bool) -> u64 {
    if !ev.is_indel {
        return SNP_PAD;
    }
    if inject_target && is_target(ev) {
        return INDEL_PAD + JAVA_STR;
    }
    match longest_str_len_at_variant(ctx, ev.start, ev.end) {
        Some(n) => INDEL_PAD + n as u64,
        None => INDEL_PAD,
    }
}

struct Replay {
    after_target: (u64, u64),
    target_padding: u64,
    after_next: (u64, u64),
    next_padding: u64,
    final_span: (u64, u64),
    injections: usize,
}

fn replay(
    events: &[Event],
    snap: &TrimCalcSnap,
    ctx: &ReferenceContext,
    inject_target: bool,
) -> Replay {
    let mut min_start = events.iter().map(|e| e.start).min().expect("events");
    let mut max_end = events.iter().map(|e| e.end).max().expect("events");
    let mut after_target = None;
    let mut target_padding = 0;
    let mut after_next = None;
    let mut next_padding = 0;
    let mut injections = 0;
    let mut seen_target = false;
    for ev in events {
        let pad = padding(ev, ctx, inject_target);
        if is_target(ev) {
            target_padding = pad;
            if inject_target {
                injections += 1;
            }
        }
        min_start = min_start.min(ev.start.saturating_sub(pad).max(1));
        max_end = max_end.max(ev.end.saturating_add(pad));
        if is_target(ev) {
            after_target = Some((min_start, max_end));
            seen_target = true;
        }
        if seen_target && is_next_indel(ev) && after_next.is_none() {
            after_next = Some((min_start, max_end));
            next_padding = pad;
        }
    }
    let final_span = (
        snap.assembly_padded_start.max(min_start),
        snap.assembly_padded_end.min(max_end),
    );
    Replay {
        after_target: after_target.expect("A>AT"),
        target_padding,
        after_next: after_next.expect("T>TGTTTG"),
        next_padding,
        final_span,
        injections,
    }
}

#[test]
fn forensic_6r293_str_padding_counterfactual() {
    gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r287_global_hap_floor(None, None);
    gatk_haplotypecaller::hc_genotyping_engine::set_forensic_6r288_haplotype_substitute(
        None, None, 0, None,
    );
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_target(None, 0);
    gatk_haplotypecaller::likelihood_engine::set_forensic_6r289_substitute(None, None);

    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam_path = root.join(BAM_REL);
    if !ref_fasta.is_file() || !bam_path.is_file() {
        eprintln!("skip: missing chr20_tiny BAM/REF");
        return;
    }
    let dict = SequenceDictionary::from_fasta_path(&ref_fasta).expect("dict");
    let specs = parse_intervals_cli_string(&dict, INTERVAL).expect("interval");
    let walk = traverse_assembly_region_walker(
        &dict,
        &specs,
        &ref_fasta,
        &bam_path,
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
    begin_trim_calc_observe();
    let _outcome = HaplotypeCallerEngine::call_region(
        region,
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let snap = take_trim_calc_snap().expect("trim snap");
    assert_eq!(
        (snap.padded_start, snap.padded_end),
        JAVA_SPAN,
        "6R.310 production trimmed span"
    );

    let mut cache = ReferenceWindowCache::new(&ref_fasta, 4);
    let ctx = ReferenceContext::from_interval(
        &dict,
        &mut cache,
        "20",
        snap.assembly_padded_start,
        snap.assembly_padded_end,
    )
    .expect("ref");
    let events = overlapping_in_trim_order(&snap);
    assert_eq!(events.iter().filter(|e| is_target(e)).count(), 1);
    assert_eq!(events.iter().filter(|e| is_next_indel(e)).count(), 1);

    let baseline = replay(&events, &snap, &ctx, false);
    assert_eq!(baseline.injections, 0);
    assert_eq!(baseline.target_padding, INDEL_PAD);
    assert_eq!(baseline.after_target, (29_455_569, 29_455_719));
    assert_eq!(
        baseline.final_span, PRE_6R310_SPAN,
        "pre-6R.310 helper replay"
    );

    let cf = replay(&events, &snap, &ctx, true);
    assert_eq!(cf.injections, 1);
    assert_eq!(cf.target_padding, INDEL_PAD + JAVA_STR);
    assert_eq!(cf.next_padding, INDEL_PAD);

    assert_eq!(cf.after_target, JAVA_SPAN);
    assert_eq!(cf.after_next, JAVA_SPAN);
    assert_eq!(cf.final_span, JAVA_SPAN);
    let classification = "STR_PADDING_COUNTERFACTUAL_REPRODUCES_JAVA_SPAN";

    kv(
        "injection_point",
        "padding assignment for 29455644 A>AT only",
    );
    kv("original_rust_padding", "75");
    kv("counterfactual_str_contribution", "9");
    kv(
        "counterfactual_effective_padding",
        cf.target_padding.to_string(),
    );
    kv("after_29455644", span(cf.after_target));
    kv("after_29455649", span(cf.after_next));
    kv("final", span(cf.final_span));
    kv("java_target", span(JAVA_SPAN));
    kv("classification", classification);
    kv(
        "production_change",
        "6R.310 live span is the former counterfactual; this replay still uses the pre-6R.310 helper",
    );
}
