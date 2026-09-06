//! 6R.100 holdout: first divergence between post-kernel max and filter-time max_ll.
//!
//! Skipped unless `HOLDOUT_6R100=1`. Coordinate-free contract lives in
//! `forensic_6r100_filter_input_likelihood_pipeline_contract`.
//!
//! ```text
//! HOLDOUT_6R100=1 cargo test -p gatk-haplotypecaller --test holdout_6r100_filter_input_likelihood -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::{
    begin_likelihood_pipeline_observe, begin_poorly_modeled_observe, call_disposition,
    flatten_assembly_regions, take_likelihood_pipeline_cells, take_likelihood_pipeline_snaps,
    take_poorly_modeled_cells, take_poorly_modeled_observe, traverse_assembly_region_walker,
    try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs,
    HaplotypeCallerEngine, LikelihoodPipelineCell, PoorlyModeledObserveCell, ReadFilterParams,
    WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
};
use serde_json::json;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "20:29455000-29456500";
const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const POS_SNP: u64 = 29_456_344;

const READS: &[(&str, u16)] = &[
    ("HISEQ1:11:H8GV6ADXX:2:2216:2203:76921", 147),
    ("HISEQ1:13:H8G92ADXX:1:1111:12251:89078", 83),
    ("HISEQ1:9:H8962ADXX:1:1112:19265:60083", 99),
    ("HWI-D00360:5:H814YADXX:1:1202:11051:34179", 147),
    ("HWI-D00360:5:H814YADXX:1:2207:10890:76583", 147),
    ("HWI-D00360:5:H814YADXX:2:1102:2154:52493", 163),
    ("HWI-D00360:6:H81VLADXX:2:1104:15554:2818", 83),
    ("HWI-D00360:6:H81VLADXX:2:1202:18367:85709", 163),
    ("HWI-D00360:7:H88WKADXX:1:1116:9273:30844", 83),
    ("HWI-D00360:8:H88U0ADXX:1:2108:16806:75328", 163),
    ("HISEQ1:13:H8G92ADXX:1:1205:16330:83279", 163),
    ("HISEQ1:9:H8962ADXX:2:1212:17767:73796", 83),
    ("HWI-D00360:5:H814YADXX:2:2103:4936:45407", 83),
    ("HWI-D00360:6:H81VLADXX:1:1103:1948:22968", 147),
    ("HWI-D00360:6:H81VLADXX:1:1210:4156:72506", 83),
    ("HWI-D00360:7:H88WKADXX:1:2111:4466:65743", 147),
    ("HWI-D00360:7:H88WKADXX:1:2203:20480:101193", 163),
    ("HWI-D00360:7:H88WKADXX:2:1214:6938:52704", 83),
    ("HWI-D00360:8:H88U0ADXX:1:1205:11075:4786", 147),
    ("HWI-D00360:8:H88U0ADXX:1:1213:18559:65935", 163),
    ("HWI-D00360:8:H88U0ADXX:1:2213:15618:11579", 163),
    ("HWI-D00360:8:H88U0ADXX:2:1213:15376:17578", 163),
];

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn java_keep(max_ll: f64) -> bool {
    !(max_ll < -8.0)
}

fn stage_seq(cells: &[LikelihoodPipelineCell], stage: &str) -> u32 {
    cells
        .iter()
        .filter(|c| c.stage == stage)
        .map(|c| c.seq)
        .min()
        .unwrap_or(0)
}

fn by_read_hap(cells: &[LikelihoodPipelineCell], stage: &str) -> HashMap<(String, u16, u64), f64> {
    let seq = stage_seq(cells, stage);
    let mut map = HashMap::new();
    for c in cells.iter().filter(|c| c.stage == stage && c.seq == seq) {
        map.insert((c.qname.clone(), c.flags, c.hap_fnv), c.log10_likelihood);
    }
    map
}

fn max_by_read(map: &HashMap<(String, u16, u64), f64>) -> HashMap<(String, u16), (f64, u64)> {
    let mut best: HashMap<(String, u16), (f64, u64)> = HashMap::new();
    for ((q, f, h), v) in map {
        let e = best
            .entry((q.clone(), *f))
            .or_insert((f64::NEG_INFINITY, 0));
        if *v > e.0 {
            *e = (*v, *h);
        }
    }
    best
}

fn cols_by_read(map: &HashMap<(String, u16, u64), f64>) -> HashMap<(String, u16), BTreeSet<u64>> {
    let mut out: HashMap<(String, u16), BTreeSet<u64>> = HashMap::new();
    for (q, f, h) in map.keys() {
        out.entry((q.clone(), *f)).or_default().insert(*h);
    }
    out
}

fn filter_map(cells: &[PoorlyModeledObserveCell], pass: u32) -> HashMap<(String, u16, u64), f64> {
    let mut map = HashMap::new();
    for c in cells.iter().filter(|c| c.pass == pass) {
        map.insert((c.qname.clone(), c.flags, c.hap_fnv), c.log10_likelihood);
    }
    map
}

fn bit_diff(
    a: &HashMap<(String, u16, u64), f64>,
    b: &HashMap<(String, u16, u64), f64>,
    reads: &[(&str, u16)],
) -> (usize, usize, usize, f64) {
    let mut n = 0usize;
    let mut eq = 0usize;
    let mut diff = 0usize;
    let mut max_abs = 0.0f64;
    for &(q, flags) in reads {
        let ha: BTreeSet<u64> = a
            .keys()
            .filter(|(qq, ff, _)| qq == q && *ff == flags)
            .map(|k| k.2)
            .collect();
        let hb: BTreeSet<u64> = b
            .keys()
            .filter(|(qq, ff, _)| qq == q && *ff == flags)
            .map(|k| k.2)
            .collect();
        for h in ha.intersection(&hb) {
            let av = a[&(q.to_string(), flags, *h)];
            let bv = b[&(q.to_string(), flags, *h)];
            n += 1;
            if av.to_bits() == bv.to_bits() {
                eq += 1;
            } else {
                diff += 1;
                max_abs = max_abs.max((av - bv).abs());
            }
        }
    }
    (n, eq, diff, max_abs)
}

fn fingerprint(map: &HashMap<(String, u16, u64), f64>, q: &str, flags: u16) -> BTreeMap<u64, u64> {
    let mut fp = BTreeMap::new();
    for ((qq, ff, h), v) in map {
        if qq == q && *ff == flags {
            fp.insert(*h, v.to_bits());
        }
    }
    fp
}

/// Restrict `src` to the haplotype keys present in `dst`. Empty if `src` is missing any of them.
fn project_onto(src: &BTreeMap<u64, u64>, dst_keys: &BTreeSet<u64>) -> Option<BTreeMap<u64, u64>> {
    let mut out = BTreeMap::new();
    for h in dst_keys {
        let bits = src.get(h)?;
        out.insert(*h, *bits);
    }
    Some(out)
}

#[test]
fn holdout_6r100_filter_input_likelihood_pipeline() {
    if std::env::var("HOLDOUT_6R100").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R100=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    assert!(ref_fasta.is_file() && bam.is_file());
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

    begin_likelihood_pipeline_observe();
    begin_poorly_modeled_observe();
    let outcome = HaplotypeCallerEngine::call_region(
        covering[0],
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let cells = take_likelihood_pipeline_cells();
    let snaps = take_likelihood_pipeline_snaps();
    let filter_rows = take_poorly_modeled_observe();
    let filter_cells = take_poorly_modeled_cells();
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
        .expect("canonical T/C");

    for s in &snaps {
        eprintln!(
            "snap seq={} stage={} n_reads={} n_haps={} n_ll={}",
            s.seq, s.stage, s.n_reads, s.n_haps, s.n_ll_entries
        );
    }

    let post = by_read_hap(&cells, "post_kernel");
    let compact = by_read_hap(&cells, "compaction");
    let norm = by_read_hap(&cells, "normalize");
    let last_pass = filter_rows.iter().map(|r| r.pass).max().unwrap_or(0);
    let filt = filter_map(&filter_cells, last_pass);
    let filter_max: HashMap<(String, u16), f64> = filter_rows
        .iter()
        .filter(|r| r.pass == last_pass)
        .map(|r| ((r.qname.clone(), r.flags), r.max_ll))
        .collect();
    let filter_row_index: HashMap<(String, u16), usize> = filter_rows
        .iter()
        .filter(|r| r.pass == last_pass)
        .map(|r| ((r.qname.clone(), r.flags), r.row_index))
        .collect();
    let post_qname_by_index: HashMap<usize, (String, u16)> = cells
        .iter()
        .filter(|c| c.stage == "post_kernel" && c.seq == stage_seq(&cells, "post_kernel"))
        .map(|c| (c.read_index, (c.qname.clone(), c.flags)))
        .collect();

    let post_max = max_by_read(&post);
    let compact_max = max_by_read(&compact);
    let norm_max = max_by_read(&norm);
    let filt_max = max_by_read(&filt);
    let post_cols = cols_by_read(&post);
    let compact_cols = cols_by_read(&compact);
    let norm_cols = cols_by_read(&norm);
    let filt_cols = cols_by_read(&filt);

    let (cn, ceq, cdiff, cabs) = bit_diff(&post, &compact, READS);
    let (nn, neq, ndiff, nabs) = bit_diff(&compact, &norm, READS);
    let (fn_, feq, fdiff, fabs) = bit_diff(&norm, &filt, READS);

    let mut n_max_post_vs_compact = 0usize;
    let mut n_max_compact_vs_norm = 0usize;
    let mut n_max_norm_vs_filt = 0usize;
    let mut n_winner_post_vs_filt = 0usize;
    let mut n_col_post_vs_compact = 0usize;
    let mut n_col_norm_vs_filt = 0usize;
    let mut n_row_remap = 0usize;
    let mut n_keep_flip = 0usize;
    let mut n_index_qname_mismatch = 0usize;
    let mut max_transform = 0.0f64;
    let mut rows = Vec::new();

    eprintln!(
        "QNAME\tfilt_idx\tscored_at_idx\tpost_max\tcomp_max\tnorm_max\tobs_max\tfilt_cell_max\tpost_win\tfilt_win\tpost_ncols\tfilt_ncols\trow_src\tdelta"
    );

    for &(q, flags) in READS {
        let k = (q.to_string(), flags);
        let (pmax, pwin) = post_max.get(&k).copied().unwrap_or((f64::NEG_INFINITY, 0));
        let (cmax, _cwin) = compact_max
            .get(&k)
            .copied()
            .unwrap_or((f64::NEG_INFINITY, 0));
        let (nmax, _nwin) = norm_max.get(&k).copied().unwrap_or((f64::NEG_INFINITY, 0));
        let (fmax, fwin) = filt_max.get(&k).copied().unwrap_or((f64::NEG_INFINITY, 0));
        let obs = filter_max.get(&k).copied().unwrap_or(f64::NEG_INFINITY);
        let pc = post_cols.get(&k).map(|s| s.len()).unwrap_or(0);
        let cc = compact_cols.get(&k).map(|s| s.len()).unwrap_or(0);
        let nc = norm_cols.get(&k).map(|s| s.len()).unwrap_or(0);
        let fc = filt_cols.get(&k).map(|s| s.len()).unwrap_or(0);
        if pmax.to_bits() != cmax.to_bits() {
            n_max_post_vs_compact += 1;
        }
        if cmax.to_bits() != nmax.to_bits() {
            n_max_compact_vs_norm += 1;
        }
        if nmax.to_bits() != fmax.to_bits() || (nmax - obs).abs() > 1e-12 {
            n_max_norm_vs_filt += 1;
        }
        if pwin != fwin {
            n_winner_post_vs_filt += 1;
        }
        if pc != cc {
            n_col_post_vs_compact += 1;
        }
        if nc != fc {
            n_col_norm_vs_filt += 1;
        }
        let delta = obs - pmax;
        max_transform = max_transform.max(delta.abs());
        if java_keep(pmax) != java_keep(obs) {
            n_keep_flip += 1;
        }
        let filt_idx = filter_row_index.get(&k).copied();
        let scored_at_idx = filt_idx.and_then(|i| post_qname_by_index.get(&i).cloned());
        let index_mismatch = match (&scored_at_idx, filt_idx) {
            (Some((sq, sf)), Some(_)) => sq != q || *sf != flags,
            _ => false,
        };
        if index_mismatch {
            n_index_qname_mismatch += 1;
        }
        let scored_label = scored_at_idx
            .as_ref()
            .map(|(sq, sf)| format!("{sq}:{sf}"))
            .unwrap_or_else(|| "-".to_string());
        let filt_idx_label = filt_idx
            .map(|i| i.to_string())
            .unwrap_or_else(|| "-".to_string());

        let fp_filt = fingerprint(&filt, q, flags);
        let filt_keys: BTreeSet<u64> = fp_filt.keys().copied().collect();
        let mut src = "self".to_string();
        if !fp_filt.is_empty() {
            let fp_self = fingerprint(&post, q, flags);
            let self_proj = project_onto(&fp_self, &filt_keys);
            if self_proj.as_ref() != Some(&fp_filt) {
                let mut found = false;
                for &(oq, of) in READS {
                    if oq == q && of == flags {
                        continue;
                    }
                    let fp_o = fingerprint(&post, oq, of);
                    if project_onto(&fp_o, &filt_keys).as_ref() == Some(&fp_filt) {
                        src = format!("{oq}:{of}");
                        n_row_remap += 1;
                        found = true;
                        break;
                    }
                }
                if !found {
                    src = "unmatched".to_string();
                    n_row_remap += 1;
                }
            }
        }

        eprintln!(
            "{q}\t{filt_idx_label}\t{scored_label}\t{pmax}\t{cmax}\t{nmax}\t{obs}\t{fmax}\t{pwin:x}\t{fwin:x}\t{pc}/{cc}/{nc}\t{fc}\t{src}\t{delta}"
        );
        rows.push(json!({
            "qname": q,
            "flags": flags,
            "filter_row_index": filt_idx,
            "scored_qname_at_filter_index": scored_label,
            "index_qname_mismatch": index_mismatch,
            "post_max": pmax,
            "compact_max": cmax,
            "norm_max": nmax,
            "filter_obs_max": obs,
            "filter_cell_max": fmax,
            "post_win": format!("{:x}", pwin),
            "filter_win": format!("{:x}", fwin),
            "post_ncols": pc,
            "compact_ncols": cc,
            "norm_ncols": nc,
            "filter_ncols": fc,
            "delta_filter_minus_post": delta,
            "row_source": src,
            "post_keep": java_keep(pmax),
            "filter_keep": java_keep(obs),
        }));
    }

    let compact_col_only = cdiff == 0 && n_col_post_vs_compact > 0;
    let classification = if n_index_qname_mismatch > 0 {
        "LIKELIHOOD_OBJECT_LIFECYCLE"
    } else if cdiff > 0 && !compact_col_only {
        "LIKELIHOOD_VALUE_TRANSFORMATION"
    } else if ndiff > 0 && n_max_compact_vs_norm > 0 {
        "LIKELIHOOD_VALUE_TRANSFORMATION"
    } else {
        "LIKELIHOOD_OBJECT_LIFECYCLE"
    };
    let first_boundary = if n_index_qname_mismatch > 0 {
        "filter_evidence_list"
    } else if cdiff > 0 || n_max_post_vs_compact > 0 {
        "compaction"
    } else if n_max_compact_vs_norm > 0 {
        "normalize"
    } else {
        "filter_evidence_aligned"
    };

    eprintln!(
        "{}",
        json!({
            "classification": classification,
            "first_boundary": first_boundary,
            "snaps": snaps.iter().map(|s| json!({
                "seq": s.seq, "stage": s.stage, "n_reads": s.n_reads,
                "n_haps": s.n_haps, "n_ll": s.n_ll_entries
            })).collect::<Vec<_>>(),
            "post_vs_compact": {"n": cn, "eq": ceq, "diff": cdiff, "max_abs": cabs},
            "compact_vs_norm": {"n": nn, "eq": neq, "diff": ndiff, "max_abs": nabs},
            "norm_vs_filter": {"n": fn_, "eq": feq, "diff": fdiff, "max_abs": fabs},
            "n_max_post_vs_compact": n_max_post_vs_compact,
            "n_max_compact_vs_norm": n_max_compact_vs_norm,
            "n_max_norm_vs_filter": n_max_norm_vs_filt,
            "n_winner_post_vs_filter": n_winner_post_vs_filt,
            "n_col_post_vs_compact": n_col_post_vs_compact,
            "n_col_norm_vs_filter": n_col_norm_vs_filt,
            "n_row_remap": n_row_remap,
            "n_index_qname_mismatch": n_index_qname_mismatch,
            "n_keep_flip": n_keep_flip,
            "max_transform_abs": max_transform,
            "vcf_ad": vcf.samples.first().map(|s| s.ad.clone()),
            "rows": rows,
        })
    );

    assert_eq!(READS.len(), 22);
    // Compaction copies retained cells (6R.96). Floor-normalize does not change max.
    assert_eq!(cdiff, 0);
    assert_eq!(n_max_post_vs_compact, 0);
    assert_eq!(n_max_compact_vs_norm, 0);
    // Java contract: filter evidence i is the record scored as read_index i.
    assert_eq!(n_index_qname_mismatch, 0);
    assert_eq!(n_keep_flip, 0);
    for r in &rows {
        let post = r["post_max"].as_f64().unwrap();
        let filt = r["filter_obs_max"].as_f64().unwrap();
        assert_eq!(
            post.to_bits(),
            filt.to_bits(),
            "{} filter max must equal post-kernel max",
            r["qname"]
        );
    }
    let _ = (
        classification,
        first_boundary,
        compact_col_only,
        n_row_remap,
        n_winner_post_vs_filt,
        n_col_norm_vs_filt,
        outcome,
    );
}
