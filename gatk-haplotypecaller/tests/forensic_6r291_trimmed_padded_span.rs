//! 6R.291: first divergence inside the trimmed-span calculation.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r291_trimmed_padded_span -- --nocapture --test-threads=1
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    begin_trim_calc_observe, call_disposition, flatten_assembly_regions, take_trim_calc_snap,
    traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_455_649;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R291\t{key}\t{}", value.as_ref());
}

#[test]
fn forensic_6r291_trimmed_padded_span() {
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
    kv(
        "active",
        format!("{}-{}", snap.active_start, snap.active_end),
    );
    kv(
        "assembly_padded",
        format!(
            "{}-{}",
            snap.assembly_padded_start, snap.assembly_padded_end
        ),
    );
    kv("n_events", snap.events.len().to_string());
    kv("n_trim_vars", snap.trim_vars.len().to_string());
    let mut overlap: Vec<_> = snap
        .events
        .iter()
        .filter(|e| e.start <= snap.active_end && e.end >= snap.active_start)
        .collect();
    overlap.sort_by_key(|e| (e.start, e.end, e.ref_al.clone(), e.alt_al.clone()));
    kv("n_overlap_events", overlap.len().to_string());
    for (i, e) in overlap.iter().enumerate() {
        kv(
            "event",
            format!(
                "i={i}\tstart={}\tend={}\tref={}\talt={}\tindel={}",
                e.start, e.end, e.ref_al, e.alt_al, e.is_indel
            ),
        );
    }
    let mut vars: Vec<_> = snap
        .trim_vars
        .iter()
        .filter(|v| v.start <= snap.active_end && v.end >= snap.active_start)
        .collect();
    vars.sort_by_key(|v| (v.start, v.end, v.is_indel));
    kv("n_overlap_trim_vars", vars.len().to_string());
    for (i, v) in vars.iter().enumerate() {
        let extra = snap
            .events
            .iter()
            .any(|e| e.start == v.start && e.end == v.end && e.is_indel == v.is_indel);
        kv(
            "trim_var",
            format!(
                "i={i}\tstart={}\tend={}\tindel={}\tin_events={extra}",
                v.start, v.end, v.is_indel
            ),
        );
    }
    kv(
        "variant_span",
        format!("{}-{}", snap.variant_start, snap.variant_end),
    );
    kv(
        "padded",
        format!("{}-{}", snap.padded_start, snap.padded_end),
    );
    let java = std::fs::read_to_string(root.join("gatk-haplotypecaller/tests/6r291_java_trim.tsv"))
        .expect("java trim tsv");
    let java_events = java_events(&java);
    assert_eq!(overlap.len(), java_events.len());
    for (i, (r, j)) in overlap.iter().zip(java_events.iter()).enumerate() {
        assert_eq!(r.start, j.start, "event {i} start");
        assert_eq!(r.end, j.end, "event {i} end");
        assert_eq!(r.alt_al, j.alt, "event {i} alt");
        assert_eq!(r.is_indel, j.indel, "event {i} indel");
        assert_eq!(r.ref_al, j.ref_al, "event {i} ref");
    }
    let j6 = &java_events[6];
    let r6 = overlap[6];
    assert_eq!(
        (
            j6.start,
            j6.end,
            j6.ref_al.as_str(),
            j6.alt.as_str(),
            j6.indel
        ),
        (29_455_644, 29_455_644, "A", "AT", true)
    );
    assert_eq!(
        (
            r6.start,
            r6.end,
            r6.ref_al.as_str(),
            r6.alt_al.as_str(),
            r6.is_indel
        ),
        (29_455_644, 29_455_644, "A", "AT", true)
    );
    assert_eq!(j6.padding, 84);
    assert_eq!((j6.left, j6.right), (29_455_560, 29_455_728));
    // Rust STR helper returns None when the event span is one base, so this indel keeps padding 75.
    let rust_pad = 75u64;
    let rust_left = r6.start - rust_pad;
    let rust_right_here = r6.end + rust_pad;
    assert_eq!(rust_left, 29_455_569);
    assert_eq!(rust_right_here, 29_455_719);
    let later = overlap
        .iter()
        .find(|e| e.start == 29_455_649 && e.is_indel)
        .expect("29455649 indel");
    let rust_right_final = rust_right_here.max(later.end + rust_pad);
    assert_eq!(rust_right_final, 29_455_724);
    assert_eq!(
        (snap.variant_start, snap.variant_end),
        (29_455_590, 29_455_703)
    );
    assert_eq!(
        (snap.active_start, snap.active_end),
        (29_455_560, 29_455_744)
    );
    assert_eq!(
        (snap.assembly_padded_start, snap.assembly_padded_end),
        (29_455_460, 29_455_844)
    );
    assert_eq!(
        (snap.padded_start, snap.padded_end),
        (rust_left, rust_right_final)
    );
    kv(
        "first_input",
        "IDENTICAL overlapping events; padding constants 20/75/75",
    );
    kv(
        "first_op",
        "indel padding at 29455644 A>AT: Java TandemRepeat longest=9 padding=84; Rust one-base span padding=75",
    );
    kv("java_after_op", "29455560-29455728");
    kv("rust_after_op", "29455569-29455719");
    kv("rust_final_right_from_next_shared_indel", "29455724");
    kv("classification", "TRIM_STR_PADDING_DIVERGENCE");
    kv("production_change", "NONE");
    kv("counterfactual", "NOT_RUN");
}

struct JavaEvent {
    start: u64,
    end: u64,
    ref_al: String,
    alt: String,
    indel: bool,
    padding: u64,
    left: u64,
    right: u64,
}

fn java_events(text: &str) -> Vec<JavaEvent> {
    let mut out = Vec::new();
    for line in text.lines() {
        let mut f = line.split('\t');
        if f.next() != Some("6R291") || f.next() != Some("event") {
            continue;
        }
        let rest: Vec<&str> = f.collect();
        let field = |k: &str| {
            rest.iter()
                .find_map(|p| p.strip_prefix(&format!("{k}=")))
                .unwrap()
                .to_string()
        };
        out.push(JavaEvent {
            start: field("start").parse().unwrap(),
            end: field("end").parse().unwrap(),
            ref_al: field("ref"),
            alt: field("alt"),
            indel: field("indel") == "true",
            padding: field("padding").parse().unwrap(),
            left: field("left").parse().unwrap(),
            right: field("right").parse().unwrap(),
        });
    }
    out
}
