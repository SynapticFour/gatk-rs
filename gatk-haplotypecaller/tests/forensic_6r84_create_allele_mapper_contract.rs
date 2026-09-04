//! 6R.84 coordinate-free: `createAlleleMapper` REF membership.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//! `AssemblyBasedCallerUtils.createAlleleMapper`:
//!
//! ```text
//! spanning = EventMap.getOverlappingEvents(loc)
//! if spanning.isEmpty()                         → REF
//! if event.start == loc and allele in mergedVC  → that ALT
//! if event.start < loc and emitSpanningDels     → SPAN_DEL (`*`)
//! if event.start < loc and !emitSpanningDels    → REF
//! unmatched overlapping events at loc           → no pool (not REF)
//! ```
//!
//! Rust colocated merge previously dumped every haplotype not claimed by an ALT
//! into REF. That leftover dump is the first mapper divergence.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r84_create_allele_mapper_contract
//! HOLDOUT_6R84=1 cargo test -p gatk-haplotypecaller --test forensic_6r84_create_allele_mapper_contract live_ -- --nocapture
//! ```

use gatk_haplotypecaller::event_map::{
    overlapping_events, remap_alt_onto_longer_ref, VariationEvent,
};
use gatk_haplotypecaller::genome_loc::GenomePosition;
use std::collections::{HashMap, HashSet};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum JavaPool {
    Ref,
    Alt(usize),
    SpanDel,
    Unassigned,
}

/// GATK 4.4 `createAlleleMapper` walk for one haplotype (first matching pool).
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

/// Pre-6R.84 Rust colocated merge: unmatched leftovers become REF.
fn rust_leftover_dump_pool(java: JavaPool) -> JavaPool {
    match java {
        JavaPool::Unassigned | JavaPool::SpanDel => JavaPool::Ref,
        other => other,
    }
}

fn ev(start: u64, r: &str, a: &str) -> VariationEvent {
    VariationEvent::from_alleles("chr", start, r, a)
}

#[test]
fn forensic_6r84_empty_overlap_is_ref() {
    let alts = vec!["T".to_string(), "CG".to_string()];
    assert_eq!(java_mapper_pool(&[], 100, "TG", &alts, true), JavaPool::Ref);
}

#[test]
fn forensic_6r84_matching_snp_and_deletion_are_alts() {
    let alts = vec!["T".to_string(), "CG".to_string()];
    assert_eq!(
        java_mapper_pool(&[ev(100, "TG", "T")], 100, "TG", &alts, true),
        JavaPool::Alt(0)
    );
    assert_eq!(
        java_mapper_pool(&[ev(100, "T", "C")], 100, "TG", &alts, true),
        JavaPool::Alt(1)
    );
}

#[test]
fn forensic_6r84_unmatched_at_loc_is_not_ref() {
    let alts = vec!["T".to_string(), "CG".to_string()];
    let java = java_mapper_pool(&[ev(100, "T", "A")], 100, "TG", &alts, true);
    assert_eq!(java, JavaPool::Unassigned);
    assert_eq!(
        rust_leftover_dump_pool(java),
        JavaPool::Ref,
        "leftover dump is the first mapper divergence"
    );
}

#[test]
fn forensic_6r84_spanning_elsewhere_with_emit_is_span_del_not_ref() {
    let alts = vec!["T".to_string(), "CG".to_string()];
    // Deletion starting before loc whose inclusive end covers loc.
    let spanning = ev(98, "AAAAA", "A");
    assert!(spanning.end_1based.get() >= 100);
    let java = java_mapper_pool(&[spanning], 100, "TG", &alts, true);
    assert_eq!(java, JavaPool::SpanDel);
    assert_eq!(rust_leftover_dump_pool(java), JavaPool::Ref);
}

#[test]
fn forensic_6r84_spanning_elsewhere_without_emit_is_ref() {
    let alts = vec!["T".to_string(), "CG".to_string()];
    let spanning = ev(98, "AAAAA", "A");
    assert_eq!(
        java_mapper_pool(&[spanning], 100, "TG", &alts, false),
        JavaPool::Ref
    );
}

#[test]
fn forensic_6r84_get_overlapping_events_matches_java_del_ins_drop() {
    let loc = 110u64;
    let del = ev(108, "AAA", "A"); // simple del ending at 110
    let ins = ev(110, "A", "AT");
    assert_eq!(del.end_1based.get(), loc);
    let both = overlapping_events(&[del.clone(), ins.clone()], loc);
    assert_eq!(
        both.len(),
        1,
        "Java drops the deletion ending at loc when an insertion also overlaps"
    );
    assert_eq!(both[0].alt_allele, "AT");
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

#[test]
fn live_create_allele_mapper_ref_pool() {
    if std::env::var("HOLDOUT_6R84").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R84=1");
        return;
    }
    use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
    use gatk_haplotypecaller::event_map::build_per_haplotype_variation_events;
    use gatk_haplotypecaller::hc_allele_mapping::{
        create_allele_mapper_with_events, hap_base_at_ref_locus,
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

    let prod_seqs: Vec<String> = outcome
        .assembly
        .haplotypes
        .iter()
        .map(|h| String::from_utf8_lossy(&h.bases).into_owned())
        .collect();
    let prod_set: HashSet<&str> = prod_seqs.iter().map(|s| s.as_str()).collect();
    let kernel_set: HashSet<&str> = java_reg.haps.iter().map(|s| s.as_str()).collect();
    let common: HashSet<&str> = kernel_set.intersection(&prod_set).copied().collect();
    let java_only: Vec<&str> = kernel_set.difference(&prod_set).copied().collect();
    let rust_only: Vec<&str> = prod_set.difference(&kernel_set).copied().collect();
    eprintln!(
        "6R.84 70→68 COMMON={} JAVA_ONLY={} RUST_ONLY={} prod={}",
        common.len(),
        java_only.len(),
        rust_only.len(),
        prod_seqs.len()
    );
    assert_eq!(prod_seqs.len(), 68);
    assert_eq!(common.len(), 68);
    assert_eq!(java_only.len(), 2);
    assert_eq!(rust_only.len(), 0);

    for s in &java_only {
        let col = java_reg.haps.iter().position(|h| h == s);
        let dup_of_kept = prod_seqs.iter().any(|p| p == *s);
        eprintln!(
            "6R.84 FILTERED_COL hash={:016x} len={} kernel_col={:?} unique_vs_68={} is_dup_seq={}",
            fnv1a64(s),
            s.len(),
            col,
            !dup_of_kept,
            dup_of_kept
        );
    }

    let haps = &outcome.assembly.haplotypes;
    let full_pad = outcome.assembly.padded_reference_start_1based();
    let full_ref = outcome.assembly.reference_bases();
    let (apply_bytes, apply_pad) = outcome.assembly.event_map_reference();
    let contig = covering[0].contig.as_str();
    let max_mnp = outcome.assembly.max_mnp_distance();
    let hap_events_full =
        build_per_haplotype_variation_events(haps, full_ref, full_pad, max_mnp, contig);
    let hap_events_apply =
        build_per_haplotype_variation_events(haps, apply_bytes, apply_pad, max_mnp, contig);

    let live = take_colocated_merge_numerics();
    let snap = live
        .iter()
        .find(|n| n.loc == POS_SNP)
        .cloned()
        .expect("colocated merge numerics");
    eprintln!(
        "6R.84 merge long_ref={} alts={:?} rust_pools={:?} java_style={:?} rust_ref_java_unassigned={} unassigned_java={} hap_sigs={:?}",
        snap.long_ref,
        snap.alts,
        snap.pool_sizes,
        snap.java_style_pool_sizes,
        snap.n_haps_rust_ref_but_java_unassigned,
        snap.n_haps_unassigned_java,
        snap.hap_event_signatures_at_loc
    );

    let long_ref = snap.long_ref.as_str();
    let alts = snap.alts.as_slice();
    let emit_spanning = true;
    let end_1based = POS_SNP.saturating_add(long_ref.len().saturating_sub(1) as u64);

    let mut rust_alt: HashSet<usize> = HashSet::new();
    for alt in alts {
        let ev = VariationEvent {
            contig: contig.to_string(),
            start_1based: GenomePosition::new_1based(POS_SNP),
            end_1based: GenomePosition::new_1based(end_1based),
            ref_allele: long_ref.to_string(),
            alt_allele: alt.clone(),
        };
        let mapping = create_allele_mapper_with_events(
            &ev,
            POS_SNP,
            haps,
            apply_pad,
            apply_bytes,
            max_mnp,
            emit_spanning,
            Some(&hap_events_full),
        );
        for hi in mapping.alt_haplotype_indices {
            rust_alt.insert(hi.get());
        }
    }
    let rust_ref: Vec<usize> = (0..haps.len()).filter(|i| !rust_alt.contains(i)).collect();

    let mut java_ref = Vec::new();
    let mut java_alt: Vec<Vec<usize>> = vec![Vec::new(); alts.len()];
    let mut java_unassigned = Vec::new();
    let mut java_span = Vec::new();
    let mut extra_ref = Vec::new();
    for i in 0..haps.len() {
        let spanning = overlapping_events(hap_events_full.events_for(i), POS_SNP);
        let spanning_apply = overlapping_events(hap_events_apply.events_for(i), POS_SNP);
        let pool = java_mapper_pool(&spanning, POS_SNP, long_ref, alts, emit_spanning);
        match pool {
            JavaPool::Ref => java_ref.push(i),
            JavaPool::Alt(ai) => java_alt[ai].push(i),
            JavaPool::Unassigned => java_unassigned.push(i),
            JavaPool::SpanDel => java_span.push(i),
        }
        let rust_leftover_ref = !rust_alt.contains(&i);
        if rust_leftover_ref && pool != JavaPool::Ref {
            extra_ref.push(i);
        }
        if pool == JavaPool::SpanDel {
            let base = hap_base_at_ref_locus(&haps[i], apply_pad, POS_SNP);
            let seq = String::from_utf8_lossy(&haps[i].bases);
            let kernel_col = java_reg.haps.iter().position(|h| h == seq.as_ref());
            eprintln!(
                "6R.84 SPAN_DEL i={} hash={:016x} len={} is_ref_flag={} kernel_col={:?} base_at_loc={:?} leftover_would_be_ref={} n_events_full={} n_events_apply={} overlap_full={:?} overlap_apply={:?} at_loc_full={:?} cigar={:?}",
                i,
                fnv1a64(seq.as_ref()),
                haps[i].bases.len(),
                haps[i].is_reference,
                kernel_col,
                base.map(|b| b as char),
                rust_leftover_ref,
                hap_events_full.events_for(i).len(),
                hap_events_apply.events_for(i).len(),
                spanning
                    .iter()
                    .map(|e| format!("{}:{}/{}", e.start_1based.get(), e.ref_allele, e.alt_allele))
                    .collect::<Vec<_>>(),
                spanning_apply
                    .iter()
                    .map(|e| format!("{}:{}/{}", e.start_1based.get(), e.ref_allele, e.alt_allele))
                    .collect::<Vec<_>>(),
                hap_events_full
                    .events_for(i)
                    .iter()
                    .filter(|e| e.start_1based.get() == POS_SNP)
                    .map(|e| format!("{}/{}", e.ref_allele, e.alt_allele))
                    .collect::<Vec<_>>(),
                haps[i]
                    .cigar
                    .as_ref()
                    .map(|c| c.to_gatk_string()),
            );
        }
        let _ = spanning_apply;
    }

    eprintln!(
        "6R.84 mapper rust_ref={} rust_alt={} java_ref={} java_alts={:?} java_unassigned={} java_span_del={} extra_ref={}",
        rust_ref.len(),
        rust_alt.len(),
        java_ref.len(),
        java_alt.iter().map(|p| p.len()).collect::<Vec<_>>(),
        java_unassigned.len(),
        java_span.len(),
        extra_ref.len()
    );

    // Aligned 136×allele matrix from Java dump using both pool definitions.
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

    let seq_to_kernel: HashMap<&str, usize> = java_reg
        .haps
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i))
        .collect();
    let rust_ref_cols: Vec<usize> = rust_ref
        .iter()
        .filter_map(|&i| seq_to_kernel.get(prod_seqs[i].as_str()).copied())
        .collect();
    let java_ref_cols: Vec<usize> = java_ref
        .iter()
        .filter_map(|&i| seq_to_kernel.get(prod_seqs[i].as_str()).copied())
        .collect();
    let mut rust_alt_cols: Vec<Vec<usize>> = vec![Vec::new(); alts.len()];
    let mut java_alt_cols = rust_alt_cols.clone();
    for (ai, pool) in java_alt.iter().enumerate() {
        java_alt_cols[ai] = pool
            .iter()
            .filter_map(|&i| seq_to_kernel.get(prod_seqs[i].as_str()).copied())
            .collect();
    }
    for i in 0..haps.len() {
        if let Some(&col) = seq_to_kernel.get(prod_seqs[i].as_str()) {
            if rust_alt.contains(&i) {
                // recover which alt via java walk first, else leftover
                let spanning = overlapping_events(hap_events_full.events_for(i), POS_SNP);
                if let JavaPool::Alt(ai) =
                    java_mapper_pool(&spanning, POS_SNP, long_ref, alts, true)
                {
                    rust_alt_cols[ai].push(col);
                }
            }
        }
    }

    fn pool_max(row: usize, h: usize, cols: &[usize], rows: &[DumpRow]) -> f64 {
        cols.iter()
            .filter_map(|&c| rows[row * h + c].lk)
            .fold(f64::NEG_INFINITY, f64::max)
    }
    let mut n_diff = 0usize;
    let mut max_abs = 0.0f64;
    let mut sum_abs = 0.0f64;
    let mut n_cells = 0usize;
    for &ri in &kept {
        let rust_lls = {
            let mut v = vec![pool_max(ri, h, &rust_ref_cols, &java_reg.rows)];
            for cols in &rust_alt_cols {
                v.push(pool_max(ri, h, cols, &java_reg.rows));
            }
            v
        };
        let java_lls = {
            let mut v = vec![pool_max(ri, h, &java_ref_cols, &java_reg.rows)];
            for cols in &java_alt_cols {
                v.push(pool_max(ri, h, cols, &java_reg.rows));
            }
            v
        };
        for (a, b) in rust_lls.iter().zip(java_lls.iter()) {
            if a.is_finite() && b.is_finite() {
                let d = (a - b).abs();
                n_cells += 1;
                sum_abs += d;
                if d > max_abs {
                    max_abs = d;
                }
                if d > 1e-4 {
                    n_diff += 1;
                }
            }
        }
    }
    let mean_abs = if n_cells == 0 {
        0.0
    } else {
        sum_abs / n_cells as f64
    };
    eprintln!(
        "6R.84 matrix 136x{} leftover_vs_java_pools differing_gt_1e4={} max_abs={:.6e} mean_abs={:.6e} cells={}",
        1 + alts.len(),
        n_diff,
        max_abs,
        mean_abs,
        n_cells
    );

    if let Some(v) = vcf {
        eprintln!(
            "6R.84 rust_vcf {}:{} PL={:?} AD={:?} QUAL={:?}",
            v.reference,
            v.alternate.join(","),
            v.samples.first().map(|s| &s.pl),
            v.samples.first().map(|s| &s.ad),
            v.quality
        );
    }
    eprintln!(
        "6R.84 AD merged={:?} remarg={:?} (PL not investigated)",
        snap.merged_ad, snap.subset_ad_remarginalized
    );

    assert_eq!(
        java_span.len(),
        6,
        "six haplotypes are Java SPAN_DEL, not unmatched-at-loc"
    );
    assert_eq!(snap.n_haps_rust_ref_but_java_unassigned, 0);
    assert_eq!(
        snap.pool_sizes, snap.java_style_pool_sizes,
        "production exclusive pools must match Java createAlleleMapper (SPAN_DEL not REF)"
    );
    assert_eq!(snap.pool_sizes[0], 35);
    assert_eq!(snap.pool_sizes[1], 6);
    assert_eq!(snap.pool_sizes[2], 21);
    if snap.pool_sizes.len() > 3 {
        assert_eq!(
            snap.pool_sizes[3], 6,
            "SPAN_DEL pool is the six 76M1D170M haplotypes"
        );
    }
    // 6R.84 leftover dump: without * in the merged list those six would have been REF.
    // After 6R.85 they occupy the * pool instead (`leftover_would_be_ref=false` above).
    assert!(
        extra_ref.is_empty() || extra_ref == java_span,
        "SPAN_DEL haplotypes must not remain as extra REF; extra_ref={extra_ref:?} java_span={java_span:?}"
    );
}
