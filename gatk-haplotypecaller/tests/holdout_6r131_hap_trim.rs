//! 6R.131 forensic dump: untrimmed vs post-trim unique haplotype sets at `20:29456196`.
//! Skipped unless `HOLDOUT_6R131=1`. Production change: NONE.
//!
//! ```text
//! HOLDOUT_6R131=1 cargo test -p gatk-haplotypecaller --test holdout_6r131_hap_trim -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    begin_hap_list_observe, begin_trim_io_observe, call_disposition, flatten_assembly_regions,
    take_hap_list_snaps, take_hap_list_trim_span, take_trim_io, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams,
    WalkerTraversalConfig,
};
use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 29_456_196;
const WORST_CELL_HAP: &str = "2d97ce18dff15697";

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn kv(key: &str, value: impl AsRef<str>) {
    println!("6R131\t{key}\t{}", value.as_ref());
}

fn unique_hashes(cols: &[gatk_haplotypecaller::HapListColumn]) -> BTreeSet<String> {
    cols.iter().map(|c| format!("{:016x}", c.fnv1a)).collect()
}

#[test]
fn holdout_6r131_hap_trim() {
    if std::env::var("HOLDOUT_6R131").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R131=1");
        return;
    }
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
    let covering = regions
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
        "active",
        format!(
            "{}:{}-{}",
            covering.contig,
            covering.start.get(),
            covering.end.get()
        ),
    );
    let args = CallRegionArgs::strict_java();
    begin_hap_list_observe();
    begin_trim_io_observe();
    let _ = HaplotypeCallerEngine::call_region(covering, &dict, &ref_fasta, &args)
        .expect("call")
        .expect("outcome");
    if let Some(span) = take_hap_list_trim_span() {
        kv(
            "trim_span",
            format!(
                "active={}-{} trim={}-{}",
                span.active_start, span.active_end, span.trim_start, span.trim_end
            ),
        );
    }
    let snaps = take_hap_list_snaps();
    kv("n_snaps", snaps.len().to_string());
    let mut by_stage: BTreeMap<&str, BTreeSet<String>> = BTreeMap::new();
    for s in &snaps {
        let uniq = unique_hashes(&s.columns);
        let n_ref = s.columns.iter().filter(|c| c.is_reference).count();
        let mut dup: BTreeMap<String, usize> = BTreeMap::new();
        for c in &s.columns {
            *dup.entry(format!("{:016x}", c.fnv1a)).or_insert(0) += 1;
        }
        let n_dup_seqs = dup.values().filter(|n| **n > 1).count();
        kv(
            "snap",
            format!(
                "stage={}\tn_cols={}\tunique={}\tn_ref={}\tdup_seqs={}",
                s.stage,
                s.n,
                uniq.len(),
                n_ref,
                n_dup_seqs
            ),
        );
        for c in &s.columns {
            kv(
                "hap",
                format!(
                    "stage={}\tidx={}\thash={:016x}\tisRef={}\tlen={}\tloc={}-{}\tcigar={}",
                    s.stage,
                    c.index,
                    c.fnv1a,
                    c.is_reference,
                    c.len,
                    c.loc_start,
                    c.loc_end,
                    c.cigar
                ),
            );
        }
        by_stage.insert(s.stage, uniq);
    }
    let io = take_trim_io();
    let mut n_dropped = 0usize;
    let mut n_kept = 0usize;
    let mut n_collapsed = 0usize;
    let mut identity = 0usize;
    for row in &io {
        if row.outcome == "dropped" {
            n_dropped += 1;
        } else if row.outcome == "kept" {
            n_kept += 1;
            if row.input_fnv == row.output_fnv {
                identity += 1;
            }
        } else if row.outcome == "collapsed" {
            n_collapsed += 1;
        }
        kv(
            "map",
            format!(
                "idx={}\tin={:016x}\tisRef={}\tlen={}\toutcome={}\tout={:016x}\toutLen={}",
                row.input_idx,
                row.input_fnv,
                row.input_is_ref,
                row.input_len,
                row.outcome,
                row.output_fnv,
                row.output_len
            ),
        );
    }
    kv(
        "trim_io_summary",
        format!(
            "n_input={}\tkept={}\tcollapsed={}\tdropped={}\tkept_identity_hash={}",
            io.len(),
            n_kept,
            n_collapsed,
            n_dropped,
            identity
        ),
    );
    if let (Some(before), Some(after)) =
        (by_stage.get("before_trim"), by_stage.get("after_trim_to"))
    {
        let common: BTreeSet<_> = before.intersection(after).cloned().collect();
        let lost: BTreeSet<_> = before.difference(after).cloned().collect();
        let gained: BTreeSet<_> = after.difference(before).cloned().collect();
        kv(
            "before_vs_after_unique",
            format!(
                "before={}\tafter={}\tcommon={}\tlost={}\tgained={}",
                before.len(),
                after.len(),
                common.len(),
                lost.len(),
                gained.len()
            ),
        );
        for h in &lost {
            kv("lost_across_trim", h.as_str());
        }
        for h in &gained {
            kv("gained_across_trim", h.as_str());
        }
        kv(
            "worst_cell_before_trim",
            before.contains(WORST_CELL_HAP).to_string(),
        );
        kv(
            "worst_cell_after_trim",
            after.contains(WORST_CELL_HAP).to_string(),
        );
    }
    if let (Some(assemble), Some(before)) =
        (by_stage.get("after_assemble"), by_stage.get("before_trim"))
    {
        kv(
            "assemble_vs_before_trim",
            format!(
                "assemble_unique={}\tbefore_unique={}\tsame={}",
                assemble.len(),
                before.len(),
                assemble == before
            ),
        );
    }
}
