//! 6R.85 coordinate-free: SPAN_DEL (`*`) as an internal genotyping allele.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//!
//! ```text
//! replaceSpanDels          → start < loc becomes (single-base REF, *)
//! makeMergedVariantContext → simpleMerge keeps * (createAlleleMapping does not extend it)
//! createAlleleMapper       → spanning.start < loc && emitSpanningDels → *
//! readLikelihoods.marginalize(alleleMapper)  → * is a likelihood column
//! calculateGLsForThisEvent → uses mergedVC alleles including *
//! calculateOutputAlleleSubset → drops unused alts (often * and the unused indel)
//! reverseTrimAlleles       → remaining SNP may emit without *
//! DepthPerAlleleBySample   → remarginalize to remaining call alleles
//! ```
//!
//! `*` is a genotyping category. It is not required in the emitted VCF.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r85_spanning_deletion_mapper_contract
//! HOLDOUT_6R85=1 cargo test -p gatk-haplotypecaller --test forensic_6r85_spanning_deletion_mapper_contract live_ -- --nocapture
//! ```

use gatk_haplotypecaller::event_map::{
    merged_alleles_for_genotyping, overlapping_events, remap_alt_onto_longer_ref, VariationEvent,
};
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::hc_allele_mapping::{replace_span_del_events, SPAN_DEL_ALLELE};
use gatk_haplotypecaller::{reverse_trim_alleles, subset_unused_alts_after_merged_genotyping};
use std::collections::{HashMap, HashSet};

fn ev(start: u64, r: &str, a: &str) -> VariationEvent {
    VariationEvent::from_alleles("chr", start, r, a)
}

/// Java `simpleMerge` allele accumulation after `replaceSpanDels` / `createAlleleMapping`.
fn java_simple_merge_alleles(events_at_loc: &[VariationEvent]) -> Vec<String> {
    let long_ref = events_at_loc
        .iter()
        .map(|e| e.ref_allele.as_str())
        .max_by_key(|r| r.len())
        .unwrap_or("")
        .to_string();
    let mut alleles = Vec::new();
    let mut seen = HashSet::new();
    seen.insert(long_ref.clone());
    alleles.push(long_ref.clone());
    for e in events_at_loc {
        if e.ref_allele == long_ref {
            if e.alt_allele != long_ref && seen.insert(e.alt_allele.clone()) {
                alleles.push(e.alt_allele.clone());
            }
        } else if e.alt_allele == SPAN_DEL_ALLELE {
            if seen.insert(SPAN_DEL_ALLELE.to_string()) {
                alleles.push(SPAN_DEL_ALLELE.to_string());
            }
        } else if let Some(remapped) =
            remap_alt_onto_longer_ref(&e.ref_allele, &e.alt_allele, &long_ref)
        {
            if remapped != long_ref && seen.insert(remapped.clone()) {
                alleles.push(remapped);
            }
        }
    }
    alleles
}

#[test]
fn forensic_6r85_replace_span_dels_converts_start_before_loc() {
    let loc = 100u64;
    let spanning = ev(98, "AAAAA", "A");
    let snp = ev(100, "T", "C");
    let replaced = replace_span_del_events(&[spanning, snp], loc, 90, b"NNNNNNNNNNTGNNNN");
    assert_eq!(replaced[0].start_1based.get(), loc);
    assert_eq!(replaced[0].alt_allele, SPAN_DEL_ALLELE);
    assert_eq!(replaced[0].ref_allele.len(), 1);
    assert_eq!(replaced[1].alt_allele, "C");
}

#[test]
fn forensic_6r85_star_is_not_extended_onto_longer_ref() {
    assert_eq!(
        remap_alt_onto_longer_ref("T", "*", "TG").as_deref(),
        Some("*")
    );
}

#[test]
fn forensic_6r85_simple_merge_retains_star_as_genotyping_allele() {
    let snp = ev(100, "T", "C");
    let del = ev(100, "TG", "T");
    let star = ev(100, "T", "*");
    let java = java_simple_merge_alleles(&[star.clone(), snp.clone(), del.clone()]);
    assert_eq!(java[0], "TG");
    assert!(
        java.iter().any(|a| a == SPAN_DEL_ALLELE),
        "Java simpleMerge keeps *: {java:?}"
    );
    assert!(java.iter().any(|a| a == "T"));
    assert!(java.iter().any(|a| a == "CG"));

    let (long_ref, alts) =
        merged_alleles_for_genotyping(&[snp, del, star], 100).expect("merged site");
    assert_eq!(long_ref, "TG");
    assert!(
        alts.iter().any(|a| a == SPAN_DEL_ALLELE),
        "Rust merged genotyping alleles must retain *: {alts:?}"
    );
    assert!(alts.iter().any(|a| a == "T"));
    assert!(alts.iter().any(|a| a == "CG"));
}

#[test]
fn forensic_6r85_no_star_event_stays_three_way() {
    let snp = ev(100, "T", "C");
    let del = ev(100, "TG", "T");
    let (long_ref, alts) = merged_alleles_for_genotyping(&[snp, del], 100).expect("merged");
    assert_eq!(long_ref, "TG");
    assert_eq!(alts, vec!["T".to_string(), "CG".to_string()]);
    assert!(!alts.iter().any(|a| a == SPAN_DEL_ALLELE));
}

#[test]
fn forensic_6r85_unused_alt_subset_drops_star_when_not_in_gt() {
    let alts = vec!["T".to_string(), "CG".to_string(), "*".to_string()];
    // 4-allele diploid: 10 GLs. SNP het is REF/CG = alleles 0/2.
    let gls = vec![0.0; 10];
    let ad = vec![10, 1, 8, 4];
    let subset =
        subset_unused_alts_after_merged_genotyping(&alts, &[0, 2], &gls, &ad).expect("subset");
    assert_eq!(subset.alt_alleles, vec!["CG".to_string()]);
    assert!(
        !subset.alt_alleles.iter().any(|a| a == SPAN_DEL_ALLELE),
        "* must not be emitted when unused"
    );
    let (trim_ref, trim_alts) = reverse_trim_alleles("TG", &subset.alt_alleles);
    assert_eq!(trim_ref, "T");
    assert_eq!(trim_alts, vec!["C".to_string()]);
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum JavaPool {
    Ref,
    Alt(usize),
    SpanDel,
    Unassigned,
}

fn java_mapper_pool(
    spanning: &[VariationEvent],
    loc: u64,
    long_ref: &str,
    alts: &[String],
    emit_spanning: bool,
) -> JavaPool {
    let loc_pos = GenomePosition::new_1based(loc);
    if spanning.is_empty() {
        return JavaPool::Ref;
    }
    let mut hit_alt: Option<usize> = None;
    for ev in spanning {
        if ev.start_1based == loc_pos {
            if ev.ref_allele.len() == long_ref.len() {
                if let Some(ai) = alts.iter().position(|a| a == &ev.alt_allele) {
                    hit_alt = Some(ai);
                }
            } else if ev.ref_allele.len() < long_ref.len() {
                if let Some(remapped) =
                    remap_alt_onto_longer_ref(&ev.ref_allele, &ev.alt_allele, long_ref)
                {
                    if let Some(ai) = alts.iter().position(|a| a == &remapped) {
                        hit_alt = Some(ai);
                    }
                }
            }
        } else if emit_spanning {
            return JavaPool::SpanDel;
        } else {
            return JavaPool::Ref;
        }
    }
    match hit_alt {
        Some(ai) => JavaPool::Alt(ai),
        None => JavaPool::Unassigned,
    }
}

#[test]
fn forensic_6r85_spanning_elsewhere_maps_to_star_pool_not_ref() {
    let alts = vec!["T".to_string(), "CG".to_string(), "*".to_string()];
    let spanning = ev(98, "AAAAA", "A");
    let java = java_mapper_pool(&[spanning], 100, "TG", &alts, true);
    assert_eq!(java, JavaPool::SpanDel);
    let star_idx = alts.iter().position(|a| a == SPAN_DEL_ALLELE).unwrap();
    assert_eq!(star_idx, 2);
}

#[derive(Clone, Debug, Eq, PartialEq, Hash)]
struct ReadFp {
    bases: String,
    bq: String,
    iq: String,
    dq: String,
    gcp: String,
}

#[derive(Clone, Debug)]
struct DumpRow {
    hap: String,
    read: ReadFp,
    lk: Option<f64>,
}

struct DumpRegion {
    haps: Vec<String>,
    reads: Vec<ReadFp>,
    rows: Vec<DumpRow>,
}

fn parse_dump(text: &str) -> Vec<DumpRegion> {
    let mut raw = Vec::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split(' ').collect();
        if parts.len() < 6 {
            continue;
        }
        let lk = parts.get(6).and_then(|s| s.parse::<f64>().ok());
        raw.push(DumpRow {
            hap: parts[0].to_string(),
            read: ReadFp {
                bases: parts[1].to_string(),
                bq: parts[2].to_string(),
                iq: parts[3].to_string(),
                dq: parts[4].to_string(),
                gcp: parts[5].to_string(),
            },
            lk,
        });
    }
    let mut regions = Vec::new();
    let mut i = 0;
    while i < raw.len() {
        let rid = &raw[i].read;
        let mut haps = Vec::new();
        let mut j = i;
        while j < raw.len() && raw[j].read == *rid {
            haps.push(raw[j].hap.clone());
            j += 1;
        }
        let h = haps.len();
        let mut k = j;
        while k + h <= raw.len() {
            let block: Vec<String> = (0..h).map(|t| raw[k + t].hap.clone()).collect();
            if block != haps {
                break;
            }
            k += h;
        }
        let mut reads = Vec::new();
        let mut t = i;
        while t < k {
            reads.push(raw[t].read.clone());
            t += h;
        }
        regions.push(DumpRegion {
            haps,
            reads,
            rows: raw[i..k].to_vec(),
        });
        i = k;
    }
    regions
}

fn region_with_motif<'a>(regions: &'a [DumpRegion], motif: &str) -> Option<&'a DumpRegion> {
    regions
        .iter()
        .find(|r| r.haps.iter().any(|h| h.contains(motif)))
}

fn fnv1a64(s: &str) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for b in s.as_bytes() {
        h ^= u64::from(*b);
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

fn java_keep_read(max_ll: f64, qualified_len: usize) -> bool {
    let max_errors = (qualified_len as f64 * 0.02).ceil().min(2.0);
    max_ll >= max_errors * -4.0
}

fn pool_max(row: usize, h: usize, cols: &[usize], rows: &[DumpRow]) -> f64 {
    cols.iter()
        .filter_map(|&c| rows[row * h + c].lk)
        .fold(f64::NEG_INFINITY, f64::max)
}

fn col_stats(java: &[f64], rust: &[f64]) -> (usize, f64, f64) {
    let mut n_diff = 0usize;
    let mut max_abs = 0.0f64;
    let mut sum_abs = 0.0f64;
    for (&j, &r) in java.iter().zip(rust.iter()) {
        let d = (j - r).abs();
        sum_abs += d;
        if d > max_abs {
            max_abs = d;
        }
        if d > 1e-4 {
            n_diff += 1;
        }
    }
    let mean = if java.is_empty() {
        0.0
    } else {
        sum_abs / java.len() as f64
    };
    (n_diff, max_abs, mean)
}

fn informative_vote(lls: &[f64]) -> Option<usize> {
    const THR: f64 = 0.2;
    let mut best_i = 0usize;
    let mut best = f64::NEG_INFINITY;
    let mut second = f64::NEG_INFINITY;
    for (i, &ll) in lls.iter().enumerate() {
        if ll > best {
            second = best;
            best = ll;
            best_i = i;
        } else if ll > second {
            second = ll;
        }
    }
    if best.is_finite() && (best - second).abs() > THR {
        Some(best_i)
    } else {
        None
    }
}

#[test]
fn live_spanning_deletion_internal_allele() {
    if std::env::var("HOLDOUT_6R85").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R85=1");
        return;
    }
    use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
    use gatk_haplotypecaller::event_map::{
        build_per_haplotype_variation_events, variation_events_at_position_from_cache,
    };
    use gatk_haplotypecaller::{
        call_disposition, flatten_assembly_regions, take_colocated_merge_numerics,
        traverse_assembly_region_walker, try_emit_call_region_variants,
        AssemblyRegionCallDisposition, CallRegionArgs, HaplotypeCallerEngine, ReadFilterParams,
        WalkerTraversalConfig, DEFAULT_STAND_EMIT_CONFIDENCE,
    };
    use std::path::Path;

    const INTERVAL: &str = "20:29455000-29456500";
    const BAM_REL: &str = "parity/giab/runs/local-pairhmm-diff/HG001.20-29455000-29456500.bam";
    const REF_REL: &str = "parity/realworld/assets/hs37d5.simple.fa";
    const JAVA_DUMP_REL: &str = "parity/giab/runs/local-pairhmm-diff/6r75_java_pairhmm_inputs.txt";
    const POS_SNP: u64 = 29_456_344;
    const MOTIF: &str = "GTGGCTCACGTCTGTAAT";

    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
    let ref_fasta = root.join(REF_REL);
    let bam = root.join(BAM_REL);
    let java_dump_path = root.join(JAVA_DUMP_REL);
    if !ref_fasta.is_file() || !bam.is_file() || !java_dump_path.is_file() {
        eprintln!("skip: live BAM/ref/java dump missing");
        return;
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
    let outcome = HaplotypeCallerEngine::call_region(
        covering[0],
        &dict,
        &ref_fasta,
        &CallRegionArgs::strict_java(),
    )
    .expect("call")
    .expect("Some");
    let emitted = try_emit_call_region_variants(
        covering[0],
        &outcome,
        "SAMPLE",
        DEFAULT_STAND_EMIT_CONFIDENCE,
    )
    .unwrap_or_default();
    let vcf = emitted.iter().find(|r| {
        r.position == POS_SNP && r.reference == "T" && r.alternate.iter().any(|a| a == "C")
    });

    let java_text = std::fs::read_to_string(&java_dump_path).expect("java dump");
    let java_regions = parse_dump(&java_text);
    let java_reg = region_with_motif(&java_regions, MOTIF).expect("java motif");
    assert_eq!(java_reg.haps.len(), 70);
    assert_eq!(java_reg.reads.len(), 153);

    let snap = take_colocated_merge_numerics()
        .into_iter()
        .find(|n| n.loc == POS_SNP)
        .expect("colocated merge numerics");

    let haps = &outcome.assembly.haplotypes;
    let full_pad = outcome.assembly.padded_reference_start_1based();
    let full_ref = outcome.assembly.reference_bases();
    let contig = covering[0].contig.as_str();
    let max_mnp = outcome.assembly.max_mnp_distance();
    let hap_events =
        build_per_haplotype_variation_events(haps, full_ref, full_pad, max_mnp, contig);
    let raw = variation_events_at_position_from_cache(&hap_events, POS_SNP, true);
    let replaced = replace_span_del_events(&raw, POS_SNP, full_pad, full_ref);
    let at_loc: Vec<VariationEvent> = replaced
        .into_iter()
        .filter(|e| e.start_1based.get() == POS_SNP)
        .collect();
    let java_internal = java_simple_merge_alleles(&at_loc);
    let rust_internal: Vec<String> = std::iter::once(snap.long_ref.clone())
        .chain(snap.alts.iter().cloned())
        .collect();

    eprintln!("JAVA_INTERNAL_ALLELES = {:?}", java_internal);
    eprintln!("RUST_INTERNAL_ALLELES = {:?}", rust_internal);
    let java_final: Vec<&str> = vcf
        .map(|r| {
            let mut v = vec![r.reference.as_str()];
            v.extend(r.alternate.iter().map(|s| s.as_str()));
            v
        })
        .unwrap_or_default();
    eprintln!("RUST_FINAL_EMITTED_ALLELES = {:?}", java_final);
    eprintln!("JAVA_FINAL_EMITTED_ALLELES = [T, C] (oracle VCF)");

    assert!(
        java_internal.iter().any(|a| a == SPAN_DEL_ALLELE),
        "Java retains * after replaceSpanDels/simpleMerge: {java_internal:?}"
    );
    assert!(
        rust_internal.iter().any(|a| a == SPAN_DEL_ALLELE),
        "Rust retains * internally after mapper/merge: {rust_internal:?}"
    );
    assert!(
        !java_final.iter().any(|a| *a == SPAN_DEL_ALLELE),
        "* must not be emitted in the VCF"
    );

    let long_ref = snap.long_ref.as_str();
    let alts = snap.alts.as_slice();
    let mut java_ref = Vec::new();
    let mut java_alt: Vec<Vec<usize>> = vec![Vec::new(); alts.len()];
    let mut java_span = Vec::new();
    for i in 0..haps.len() {
        let spanning = overlapping_events(hap_events.events_for(i), POS_SNP);
        match java_mapper_pool(&spanning, POS_SNP, long_ref, alts, true) {
            JavaPool::Ref => java_ref.push(i),
            JavaPool::Alt(ai) => java_alt[ai].push(i),
            JavaPool::SpanDel => {
                java_span.push(i);
                if let Some(ai) = alts.iter().position(|a| a == SPAN_DEL_ALLELE) {
                    java_alt[ai].push(i);
                }
            }
            JavaPool::Unassigned => {}
        }
    }
    eprintln!(
        "6R.85 pools rust={:?} java_style={:?} SPAN_DEL_haplotypes={} unassigned_java={}",
        snap.pool_sizes,
        snap.java_style_pool_sizes,
        java_span.len(),
        snap.n_haps_unassigned_java
    );
    assert_eq!(java_span.len(), 6, "SPAN_DEL haplotype count");
    assert_eq!(snap.pool_sizes, vec![35, 6, 21, 6]);
    assert_eq!(snap.pool_sizes, snap.java_style_pool_sizes);
    assert_eq!(snap.n_haps_unassigned_java, 0);

    let prod_seqs: Vec<String> = haps
        .iter()
        .map(|h| String::from_utf8_lossy(&h.bases).into_owned())
        .collect();
    let seq_to_kernel: HashMap<&str, usize> = java_reg
        .haps
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();

    let h = java_reg.haps.len();
    let mut per_read_max = vec![f64::NEG_INFINITY; java_reg.reads.len()];
    for (ri, _) in java_reg.reads.iter().enumerate() {
        for j in 0..h {
            if let Some(lk) = java_reg.rows[ri * h + j].lk {
                if lk > per_read_max[ri] {
                    per_read_max[ri] = lk;
                }
            }
        }
    }
    let kept: Vec<usize> = java_reg
        .reads
        .iter()
        .enumerate()
        .filter(|(i, r)| java_keep_read(per_read_max[*i], r.bq.len().max(1)))
        .map(|(i, _)| i)
        .collect();
    assert_eq!(kept.len(), 136);

    let ref_cols: Vec<usize> = java_ref
        .iter()
        .filter_map(|&i| seq_to_kernel.get(prod_seqs[i].as_str()).copied())
        .collect();
    let mut alt_cols: Vec<Vec<usize>> = vec![Vec::new(); alts.len()];
    for (ai, pool) in java_alt.iter().enumerate() {
        alt_cols[ai] = pool
            .iter()
            .filter_map(|&i| seq_to_kernel.get(prod_seqs[i].as_str()).copied())
            .collect();
    }
    let span_cols: Vec<usize> = java_span
        .iter()
        .filter_map(|&i| seq_to_kernel.get(prod_seqs[i].as_str()).copied())
        .collect();
    assert_eq!(span_cols.len(), 6);

    let star_ai = alts.iter().position(|a| a == SPAN_DEL_ALLELE);
    let t_ai = alts.iter().position(|a| a == "T");
    let cg_ai = alts.iter().position(|a| a == "CG");
    assert!(star_ai.is_some() && t_ai.is_some() && cg_ai.is_some());

    let rust_span_cols = span_cols.clone();

    let mut java_ref_ll = Vec::new();
    let mut java_t_ll = Vec::new();
    let mut java_cg_ll = Vec::new();
    let mut java_star_ll = Vec::new();
    let mut rust_ref_ll = Vec::new();
    let mut rust_t_ll = Vec::new();
    let mut rust_cg_ll = Vec::new();
    let mut rust_star_ll = Vec::new();
    let mut n_star_best = 0usize;
    let mut n_star_informative = 0usize;
    for &ri in &kept {
        let j_ref = pool_max(ri, h, &ref_cols, &java_reg.rows);
        let j_t = t_ai
            .map(|ai| pool_max(ri, h, &alt_cols[ai], &java_reg.rows))
            .unwrap_or(f64::NEG_INFINITY);
        let j_cg = cg_ai
            .map(|ai| pool_max(ri, h, &alt_cols[ai], &java_reg.rows))
            .unwrap_or(f64::NEG_INFINITY);
        let j_star = pool_max(ri, h, &span_cols, &java_reg.rows);
        java_ref_ll.push(j_ref);
        java_t_ll.push(j_t);
        java_cg_ll.push(j_cg);
        java_star_ll.push(j_star);
        rust_ref_ll.push(j_ref);
        rust_t_ll.push(j_t);
        rust_cg_ll.push(j_cg);
        rust_star_ll.push(pool_max(ri, h, &rust_span_cols, &java_reg.rows));
        let four = [j_ref, j_t, j_cg, j_star];
        let best = four.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if (best - j_star).abs() < 1e-12 {
            n_star_best += 1;
        }
        if informative_vote(&four) == Some(3) {
            n_star_informative += 1;
        }
    }

    let (d_ref, m_ref, mean_ref) = col_stats(&java_ref_ll, &rust_ref_ll);
    let (d_t, m_t, mean_t) = col_stats(&java_t_ll, &rust_t_ll);
    let (d_cg, m_cg, mean_cg) = col_stats(&java_cg_ll, &rust_cg_ll);
    let (d_star, m_star, mean_star) = col_stats(&java_star_ll, &rust_star_ll);
    eprintln!("JAVA read×allele dimension = 136 × 4 (REF / T / CG / *)");
    eprintln!(
        "RUST read×allele dimension = 136 × {}",
        snap.pool_sizes.len()
    );
    eprintln!(
        "aligned REF differing={} max_delta={:.6e} mean_delta={:.6e}",
        d_ref, m_ref, mean_ref
    );
    eprintln!(
        "aligned T differing={} max_delta={:.6e} mean_delta={:.6e}",
        d_t, m_t, mean_t
    );
    eprintln!(
        "aligned CG differing={} max_delta={:.6e} mean_delta={:.6e}",
        d_cg, m_cg, mean_cg
    );
    eprintln!(
        "aligned SPAN_DEL differing={} max_delta={:.6e} mean_delta={:.6e} n_star_best={} n_star_informative={}",
        d_star, m_star, mean_star, n_star_best, n_star_informative
    );

    let mut ad4 = [0u32; 4];
    let mut ad_remarg_ref_cg = [0u32; 2];
    let mut ad3 = [0u32; 3];
    for i in 0..kept.len() {
        let four = [java_ref_ll[i], java_t_ll[i], java_cg_ll[i], java_star_ll[i]];
        if let Some(v) = informative_vote(&four) {
            ad4[v] += 1;
        }
        let two = [java_ref_ll[i], java_cg_ll[i]];
        if let Some(v) = informative_vote(&two) {
            ad_remarg_ref_cg[v] += 1;
        }
        let three = [java_ref_ll[i], java_t_ll[i], java_cg_ll[i]];
        if let Some(v) = informative_vote(&three) {
            ad3[v] += 1;
        }
    }
    eprintln!(
        "AD 4-way informative [REF,T,CG,*]={:?} remarg{{REF,CG}}={:?} 3-way[REF,T,CG]={:?} rust_merged={:?} rust_remarg={:?} rust_vcf={:?}",
        ad4,
        ad_remarg_ref_cg,
        ad3,
        snap.merged_ad,
        snap.subset_ad_remarginalized,
        vcf.and_then(|r| r.samples.first().and_then(|s| s.ad.clone()))
    );
    eprintln!(
        "SPAN_DEL haplotype count = 6; reads with * best = {}; reads informatively supporting SPAN_DEL = {}",
        n_star_best, n_star_informative
    );
    if let Some(v) = vcf {
        eprintln!(
            "6R.85 rust_vcf T/C GT={:?} AD={:?} PL={:?} QUAL={:?}",
            v.samples.first().map(|s| &s.gt),
            v.samples.first().map(|s| &s.ad),
            v.samples.first().map(|s| &s.pl),
            v.quality
        );
        assert!(
            !v.alternate.iter().any(|a| a == SPAN_DEL_ALLELE),
            "do not emit * merely because Java genotypes it internally"
        );
    }

    assert_eq!(d_ref, 0);
    assert_eq!(d_t, 0);
    assert_eq!(d_cg, 0);
    assert_eq!(d_star, 0);
    assert_eq!(span_cols.len(), 6);

    for i in &java_span {
        eprintln!(
            "6R.85 SPAN_DEL hap i={} hash={:016x} len={}",
            i,
            fnv1a64(prod_seqs[*i].as_str()),
            haps[*i].bases.len()
        );
    }
}
