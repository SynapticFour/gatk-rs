//! 6R.113 forensic dump: stored `A/ATG` vs haplotype CIGAR-replay `A/G` at `2:92307327`.
//! Silent unless `HOLDOUT_6R113=1`. No production algorithm change.
//!
//! ```text
//! HOLDOUT_6R113=1 P12_REFERENCE=$PWD/parity/realworld/assets/hs37d5.simple.fa \
//!   cargo test -p gatk-haplotypecaller --test holdout_6r113_eventmap_replay -- --nocapture
//! ```

use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
use gatk_haplotypecaller::cigar::{Cigar, CigarOperator};
use gatk_haplotypecaller::event_map::{
    collect_variation_events, variation_events_for_haplotype, VariationEvent,
};
use gatk_haplotypecaller::haplotype_cigar::calculate_haplotype_cigar_with_strategy;
use gatk_haplotypecaller::hc_allele_mapping::hap_base_at_ref_locus;
use gatk_haplotypecaller::{
    call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
    AssemblyRegionCallDisposition, CallRegionArgs, Haplotype, HaplotypeCallerEngine,
    ReadFilterParams, SwOverhangStrategy, SwParameters, WalkerTraversalConfig,
};
use serde_json::{json, Value};
use std::path::{Path, PathBuf};

const INTERVAL: &str = "2:92300000-92350000";
const BAM_REL: &str = "parity/realworld/na12878_20k_b37/NA12878_20k.b37.bam";
const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
const TARGET: u64 = 92_307_327;
const WIN_LO: u64 = 92_307_324;
const WIN_HI: u64 = 92_307_330;

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("..")
}

fn fnv1a64_hex(bases: &[u8]) -> String {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bases {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{h:016x}")
}

fn cigar_str(h: &Haplotype) -> String {
    h.cigar
        .as_ref()
        .map(|c| {
            c.elements
                .iter()
                .map(|e| format!("{}{}", e.length, e.operator.as_char()))
                .collect::<String>()
        })
        .unwrap_or_else(|| ".".to_string())
}

fn cigar_has_indel(h: &Haplotype) -> bool {
    h.cigar
        .as_ref()
        .is_some_and(|c| c.elements.iter().any(|e| e.operator.is_indel()))
}

fn event_key(e: &VariationEvent) -> String {
    format!("{}/{}", e.ref_allele, e.alt_allele)
}

fn events_at(events: &[VariationEvent], loc: u64) -> Vec<Value> {
    events
        .iter()
        .filter(|e| e.start_1based.get() == loc)
        .map(|e| {
            json!({
                "start": e.start_1based.get(),
                "end": e.end_1based.get(),
                "ref": e.ref_allele,
                "alt": e.alt_allele,
                "indel": e.is_indel(),
            })
        })
        .collect()
}

fn local_bases(h: &Haplotype, pad: u64, lo: u64, hi: u64) -> String {
    (lo..=hi)
        .map(|p| {
            hap_base_at_ref_locus(h, pad, p)
                .map(|b| (b as char).to_ascii_uppercase())
                .unwrap_or('-')
        })
        .collect()
}

fn ref_slice(ref_bytes: &[u8], pad: u64, lo: u64, hi: u64) -> String {
    (lo..=hi)
        .map(|p| {
            let off = p.saturating_sub(pad) as usize;
            ref_bytes
                .get(off)
                .map(|&b| (b as char).to_ascii_uppercase())
                .unwrap_or('.')
        })
        .collect()
}

/// Forced coupled-cluster CIGAR (`…M 2D 1M 2I …`) used by `force_cluster_coupled_haplotype_cigar`.
/// `match_through_anchor`: false = Match `ttc_off` (live production); true = Match `ttc_off+1`
/// (include the T of TTC so 2D starts at the deleted TC).
fn forced_2d1m2i_cigar(pad: u64, ref_len: usize, match_through_anchor: bool) -> Cigar {
    let ttc_off = WIN_LO.saturating_sub(pad) as usize;
    let lead = if match_through_anchor {
        ttc_off.saturating_add(1)
    } else {
        ttc_off
    };
    let tail = ref_len.saturating_sub(lead + 3);
    let mut c = Cigar::new();
    if lead > 0 {
        c.push(lead, CigarOperator::Match);
    }
    c.push(2, CigarOperator::Deletion);
    c.push(1, CigarOperator::Match);
    c.push(2, CigarOperator::Insertion);
    if tail > 0 {
        c.push(tail, CigarOperator::Match);
    }
    c
}

fn effective_align_start(h: &Haplotype, ref_hap: &Haplotype, pad: u64) -> usize {
    let trim_offset = ref_hap
        .genome_loc
        .map(|g| g.start_1based().saturating_sub(pad) as usize)
        .unwrap_or(0);
    let a = h.alignment_start_hap_wrt_ref;
    if trim_offset > 0 && a < trim_offset {
        trim_offset.saturating_add(a)
    } else {
        a
    }
}

fn cigar_to_string(c: &Cigar) -> String {
    c.elements
        .iter()
        .map(|e| format!("{}{}", e.length, e.operator.as_char()))
        .collect()
}

/// Walk CIGAR across genomic `lo..=hi`.
/// `ref_pos0_start` is the EventMap cursor: 0-based offset into padded `ref_bytes`.
/// Genomic pos is 1-based inclusive (`pad + ref_pos`).
fn cursor_rows(
    h: &Haplotype,
    ref_bytes: &[u8],
    pad: u64,
    ref_pos0_start: usize,
    lo: u64,
    hi: u64,
) -> Vec<Value> {
    let Some(cigar) = &h.cigar else {
        return Vec::new();
    };
    let mut ref_pos = ref_pos0_start;
    let mut hap_pos = 0usize;
    let mut rows = Vec::new();
    for el in &cigar.elements {
        match el.operator {
            CigarOperator::Match => {
                for _ in 0..el.length {
                    let genomic = pad + ref_pos as u64;
                    if genomic >= lo && genomic <= hi {
                        let rb = ref_bytes
                            .get(ref_pos)
                            .map(|&b| (b as char).to_ascii_uppercase())
                            .unwrap_or('.');
                        let hb = h
                            .bases
                            .get(hap_pos)
                            .map(|&b| (b as char).to_ascii_uppercase())
                            .unwrap_or('-');
                        rows.push(json!({
                            "genomic": genomic,
                            "ref_pos0": ref_pos,
                            "hap_pos0": hap_pos,
                            "ref": rb.to_string(),
                            "hap": hb.to_string(),
                            "op": "M",
                            "eq": rb == hb,
                        }));
                    }
                    ref_pos += 1;
                    hap_pos += 1;
                }
            }
            CigarOperator::Deletion => {
                for _ in 0..el.length {
                    let genomic = pad + ref_pos as u64;
                    if genomic >= lo && genomic <= hi {
                        let rb = ref_bytes
                            .get(ref_pos)
                            .map(|&b| (b as char).to_ascii_uppercase())
                            .unwrap_or('.');
                        rows.push(json!({
                            "genomic": genomic,
                            "ref_pos0": ref_pos,
                            "hap_pos0": hap_pos,
                            "ref": rb.to_string(),
                            "hap": "-",
                            "op": "D",
                            "eq": false,
                        }));
                    }
                    ref_pos += 1;
                }
            }
            CigarOperator::Insertion | CigarOperator::SoftClip => {
                let genomic = pad + ref_pos as u64;
                let in_win = genomic >= lo && genomic <= hi.saturating_add(1);
                if in_win {
                    let ins: String = h.bases
                        [hap_pos..hap_pos.saturating_add(el.length).min(h.bases.len())]
                        .iter()
                        .map(|&b| (b as char).to_ascii_uppercase())
                        .collect();
                    rows.push(json!({
                        "genomic_after": genomic,
                        "ref_pos0": ref_pos,
                        "hap_pos0": hap_pos,
                        "ref": "-",
                        "hap": ins,
                        "op": el.operator.as_char().to_string(),
                        "len": el.length,
                    }));
                }
                hap_pos += el.length;
            }
            CigarOperator::HardClip => {}
        }
    }
    rows
}

fn replay_keys(events: &[VariationEvent], loc: u64) -> Vec<String> {
    events
        .iter()
        .filter(|e| e.start_1based.get() == loc)
        .map(event_key)
        .collect()
}

fn replay_window(events: &[VariationEvent], lo: u64, hi: u64) -> Vec<Value> {
    events
        .iter()
        .filter(|e| e.start_1based.get() >= lo && e.start_1based.get() <= hi)
        .map(|e| {
            json!({
                "start": e.start_1based.get(),
                "end": e.end_1based.get(),
                "ref": e.ref_allele,
                "alt": e.alt_allele,
            })
        })
        .collect()
}

fn try_forced_replay(
    h: &Haplotype,
    ref_hap: &Haplotype,
    ref_bytes: &[u8],
    pad: u64,
    max_mnp: usize,
    cigar: Cigar,
) -> Value {
    let ok = cigar.read_length() == h.bases.len();
    if !ok {
        return json!({
            "cigar": cigar_to_string(&cigar),
            "read_len": cigar.read_length(),
            "hap_len": h.bases.len(),
            "ok": false,
        });
    }
    let mut forced_hap = h.clone();
    forced_hap.cigar = Some(cigar.clone());
    let replay = variation_events_for_haplotype(&forced_hap, ref_hap, ref_bytes, pad, max_mnp, "2");
    json!({
        "cigar": cigar_to_string(&cigar),
        "ok": true,
        "at_92307327": replay_keys(&replay, TARGET),
        "window_92307320_332": replay_window(&replay, 92_307_320, 92_307_332),
    })
}

#[test]
fn holdout_6r113_eventmap_replay_dump() {
    if std::env::var("HOLDOUT_6R113").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R113=1");
        return;
    }
    let root = repo_root();
    let ref_fasta = std::env::var("P12_REFERENCE")
        .map(PathBuf::from)
        .unwrap_or_else(|_| root.join(REF_REL));
    let bam = root.join(BAM_REL);
    assert!(ref_fasta.is_file(), "missing {}", ref_fasta.display());
    assert!(bam.is_file(), "missing {}", bam.display());

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
        .expect("ActiveFull covering 2:92307327");
    let args = CallRegionArgs::strict_java();
    let outcome = HaplotypeCallerEngine::call_region(covering, &dict, &ref_fasta, &args)
        .expect("call")
        .expect("outcome");
    let pad = outcome.assembly.padded_reference_start_1based();
    let ref_bytes = outcome.assembly.reference_bases();
    let haps = &outcome.assembly.haplotypes;
    let ref_hap = haps.iter().find(|h| h.is_reference).expect("ref hap");
    let max_mnp = outcome.assembly.max_mnp_distance();

    let stored = outcome.assembly.variation_events();
    let cigar_union = collect_variation_events(haps, ref_bytes, pad, "2", max_mnp);

    let stored_at = events_at(stored, TARGET);
    let stored_ttc = events_at(stored, WIN_LO);
    let cigar_union_at = events_at(&cigar_union, TARGET);
    let stored_has_atg = stored
        .iter()
        .any(|e| e.start_1based.get() == TARGET && e.ref_allele == "A" && e.alt_allele == "ATG");
    let cigar_union_has_atg = cigar_union
        .iter()
        .any(|e| e.start_1based.get() == TARGET && e.ref_allele == "A" && e.alt_allele == "ATG");
    let cigar_union_has_ag = cigar_union
        .iter()
        .any(|e| e.start_1based.get() == TARGET && e.ref_allele == "A" && e.alt_allele == "G");

    let ref_window = ref_slice(ref_bytes, pad, WIN_LO, WIN_HI);
    let sw = SwParameters::gatk_haplotype_to_reference();
    let mut hap_rows = Vec::new();
    let mut n_replay_atg = 0usize;
    let mut n_replay_ag = 0usize;
    let mut n_seq_tatg = 0usize;
    let mut substitutions = Vec::new();
    for (i, h) in haps.iter().enumerate() {
        let replay = variation_events_for_haplotype(h, ref_hap, ref_bytes, pad, max_mnp, "2");
        let at_loc = replay_keys(&replay, TARGET);
        if at_loc.iter().any(|k| k == "A/ATG") {
            n_replay_atg += 1;
        }
        if at_loc.iter().any(|k| k == "A/G") {
            n_replay_ag += 1;
        }
        let hap_window = local_bases(h, pad, WIN_LO, WIN_HI);
        let base_at = hap_base_at_ref_locus(h, pad, TARGET)
            .map(|b| (b as char).to_ascii_uppercase().to_string())
            .unwrap_or_else(|| "-".into());
        let gl = h
            .genome_loc
            .map(|g| json!([g.start_1based(), g.end_1based()]));
        let linear_window = h.genome_loc.and_then(|g| {
            let off = WIN_LO.saturating_sub(g.start_1based()) as usize;
            h.bases.get(off..off.saturating_add(7)).map(|s| {
                s.iter()
                    .map(|&b| (b as char).to_ascii_uppercase())
                    .collect::<String>()
            })
        });
        if hap_window.starts_with("TATG")
            || linear_window
                .as_deref()
                .is_some_and(|s| s.starts_with("TATG"))
        {
            n_seq_tatg += 1;
        }
        let eff = effective_align_start(h, ref_hap, pad);
        hap_rows.push(json!({
            "idx": i,
            "hash": fnv1a64_hex(&h.bases),
            "len": h.bases.len(),
            "ref_hap_len": ref_hap.bases.len(),
            "ref_bytes_len": ref_bytes.len(),
            "len_eq_ref_hap": h.bases.len() == ref_hap.bases.len(),
            "len_eq_ref_bytes": h.bases.len() == ref_bytes.len(),
            "is_reference": h.is_reference,
            "cigar": cigar_str(h),
            "cigar_has_indel": cigar_has_indel(h),
            "align_start": h.alignment_start_hap_wrt_ref,
            "effective_align_start": eff,
            "genome_loc": gl,
            "hap_window_cigar_mapped": hap_window,
            "hap_window_linear_genome_loc": linear_window,
            "base_at_92307327": base_at,
            "replay_at_loc": at_loc,
            "replay_window_92307320_332": replay_window(&replay, 92_307_320, 92_307_332),
            "cursor_eventmap_origin": cursor_rows(h, ref_bytes, pad, eff, WIN_LO, WIN_HI),
        }));

        if at_loc.iter().any(|k| k == "A/G")
            || hap_window.starts_with("TATG")
            || linear_window
                .as_deref()
                .is_some_and(|s| s.starts_with("TATG"))
        {
            let apply_len = h
                .cigar
                .as_ref()
                .map(|c| c.reference_length())
                .unwrap_or(h.bases.len());
            let apply_pad = h.genome_loc.map(|g| g.start_1based()).unwrap_or(pad);
            let cigar_79 = forced_2d1m2i_cigar(apply_pad, apply_len, false);
            let cigar_80 = forced_2d1m2i_cigar(apply_pad, apply_len, true);
            let indel_sw = calculate_haplotype_cigar_with_strategy(
                ref_bytes,
                &h.bases,
                &sw,
                SwOverhangStrategy::Indel,
            );
            let sw_ok = indel_sw
                .as_ref()
                .is_some_and(|c| c.read_length() == h.bases.len());
            let sw_keys = if sw_ok {
                let mut sw_hap = h.clone();
                sw_hap.cigar = indel_sw.clone();
                let sw_replay =
                    variation_events_for_haplotype(&sw_hap, ref_hap, ref_bytes, pad, max_mnp, "2");
                json!({
                    "at_92307327": replay_keys(&sw_replay, TARGET),
                    "window_92307320_332": replay_window(&sw_replay, 92_307_320, 92_307_332),
                })
            } else {
                json!(null)
            };
            substitutions.push(json!({
                "idx": i,
                "hash": fnv1a64_hex(&h.bases),
                "same_sequence": true,
                "live_cigar": cigar_str(h),
                "live_replay_at_92307327": at_loc,
                "live_replay_window": replay_window(&replay, 92_307_320, 92_307_332),
                "forced_apply_pad": apply_pad,
                "forced_apply_len": apply_len,
                "sub_match_ttc_off": try_forced_replay(
                    h, ref_hap, ref_bytes, pad, max_mnp, cigar_79,
                ),
                "sub_match_ttc_off_plus_1": try_forced_replay(
                    h, ref_hap, ref_bytes, pad, max_mnp, cigar_80,
                ),
                "indel_sw_cigar": indel_sw.as_ref().map(cigar_to_string),
                "indel_sw_same_as_live": indel_sw.as_ref() == h.cigar.as_ref(),
                "indel_sw_replay": sw_keys,
            }));
        }
    }

    let doc = json!({
        "holdout": "2:92307327 stored A/ATG vs replay A/G",
        "active_full": [covering.start.get(), covering.end.get()],
        "pad": pad,
        "ref_hap_len": ref_hap.bases.len(),
        "hap_count": haps.len(),
        "coordinate": {
            "genomic": "1-based inclusive VCF",
            "hap.alignment_start_hap_wrt_ref": "0-based offset into padded reference bytes",
            "EventMap Match start_pos": "ref_loc_start_1based + ref_pos + mismatch_offset (1-based genomic)",
            "EventMap Insertion start": "ref_loc_start_1based + ref_pos - 1 (preceding ref base, 1-based inclusive)",
            "rows": [
                {"pos": 92307324, "ref": ref_slice(ref_bytes, pad, 92307324, 92307324)},
                {"pos": 92307325, "ref": ref_slice(ref_bytes, pad, 92307325, 92307325)},
                {"pos": 92307326, "ref": ref_slice(ref_bytes, pad, 92307326, 92307326)},
                {"pos": 92307327, "ref": ref_slice(ref_bytes, pad, 92307327, 92307327)},
                {"pos": 92307328, "ref": ref_slice(ref_bytes, pad, 92307328, 92307328)},
                {"pos": 92307329, "ref": ref_slice(ref_bytes, pad, 92307329, 92307329)},
            ],
        },
        "object_s_stored_union": {
            "producer": "assembly.variation_events after restore_p12_cluster_genotyping_events / inject_reference_cluster_indel_events / synthesize_cluster_motif_insertions (list inject; not per-hap CIGAR harvest)",
            "at_loc": stored_at,
            "at_92307324": stored_ttc,
            "has_A_ATG": stored_has_atg,
        },
        "object_r_cigar_replay": {
            "producer": "collect_variation_events / variation_events_for_haplotype → EventMap::from_haplotype_and_reference processCigarForInitialEvents",
            "at_loc": cigar_union_at,
            "has_A_ATG": cigar_union_has_atg,
            "has_A_G": cigar_union_has_ag,
        },
        "same_event_set_stored_vs_cigar_union": stored_at == cigar_union_at,
        "ref_window_92307324_30": ref_window,
        "n_haps_replay_A_ATG": n_replay_atg,
        "n_haps_replay_A_G": n_replay_ag,
        "n_haps_window_TATG": n_seq_tatg,
        "haps": hap_rows,
        "substitution": substitutions,
    });
    println!("{}", serde_json::to_string_pretty(&doc).unwrap());

    assert!(stored_has_atg, "Object S must still list stored A/ATG");
    // 6R.114: live production CIGAR is now EventMap-compatible (80M…). Historical
    // 6R.113 counterfactual is the substitution: 79M → A/G, 80M → A/ATG.
    let mut saw_79_ag = false;
    let mut saw_80_atg = false;
    for s in &substitutions {
        if s["sub_match_ttc_off"]["at_92307327"]
            .as_array()
            .is_some_and(|a| a.iter().any(|v| v == "A/G"))
        {
            saw_79_ag = true;
        }
        if s["sub_match_ttc_off_plus_1"]["at_92307327"]
            .as_array()
            .is_some_and(|a| a.iter().any(|v| v == "A/ATG"))
        {
            saw_80_atg = true;
        }
    }
    assert!(saw_79_ag, "6R.113 historical: 79M2D1M2I replay remains A/G");
    assert!(
        saw_80_atg,
        "6R.113 historical / 6R.114 production: 80M2D1M2I replay is A/ATG"
    );
}
