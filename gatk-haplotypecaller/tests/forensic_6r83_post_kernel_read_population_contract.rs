//! 6R.83 coordinate-free: first post-PairHMM genotyping-read population.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//! `ReadLikelihoodCalculationEngine.filterPoorlyModeledEvidence` with
//! `capLikelihoods=true` (static; `--disable-dynamic-read-disqualification` default):
//!
//! ```text
//! qualifiedLen = HMM_BASE_QUALITIES length (else read length)
//! maxErrors    = min(2.0, ceil(qualifiedLen * 0.02))
//! threshold    = maxErrors * -4.0
//! drop iff maximumLikelihoodOverAllAlleles < threshold
//! ```
//!
//! Applied inside `PairHMMLikelihoodCalculationEngine.computeReadLikelihoods`
//! after `normalizeLikelihoods` and **before** optional `filterAlleles`.
//! Normalize floors losing cells; it does not change the per-read maximum.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r83_post_kernel_read_population_contract
//! HOLDOUT_6R83=1 cargo test -p gatk-haplotypecaller --test forensic_6r83_post_kernel_read_population_contract live_ -- --nocapture
//! ```

use std::collections::{HashMap, HashSet};

const EXPECTED_ERROR_RATE_PER_BASE: f64 = 0.02;
const LOG10_QUAL_PER_ERROR: f64 = -4.0;

/// Java `log10MinTrueLikelihood(expectedErrorRatePerBase, capLikelihoods=true)`.
fn java_log10_min_true_likelihood(qualified_read_len: usize) -> f64 {
    let max_errors = (qualified_read_len as f64 * EXPECTED_ERROR_RATE_PER_BASE)
        .ceil()
        .min(2.0);
    max_errors * LOG10_QUAL_PER_ERROR
}

/// Java keep predicate: drop iff `max_ll < threshold` (IEEE: keep includes equality).
fn java_keep_read(max_ll: f64, qualified_read_len: usize) -> bool {
    max_ll >= java_log10_min_true_likelihood(qualified_read_len)
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
        let mut reads = vec![raw[i].read.clone()];
        let mut k = j;
        while k + h <= raw.len() {
            let block: Vec<String> = (0..h).map(|t| raw[k + t].hap.clone()).collect();
            if block != haps {
                break;
            }
            reads.push(raw[k].read.clone());
            k += h;
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

fn per_read_max(region: &DumpRegion) -> Vec<f64> {
    let h = region.haps.len();
    region
        .reads
        .iter()
        .enumerate()
        .map(|(ri, _)| {
            let mut m = f64::NEG_INFINITY;
            for j in 0..h {
                if let Some(lk) = region.rows[ri * h + j].lk {
                    if lk > m {
                        m = lk;
                    }
                }
            }
            m
        })
        .collect()
}

fn java_kept_indices(region: &DumpRegion) -> Vec<usize> {
    let maxs = per_read_max(region);
    maxs.iter()
        .enumerate()
        .filter(|(i, &m)| java_keep_read(m, region.reads[*i].bq.len().max(1)))
        .map(|(i, _)| i)
        .collect()
}

#[test]
fn forensic_6r83_java_static_threshold_formula() {
    assert_eq!(java_log10_min_true_likelihood(1), -4.0);
    assert_eq!(java_log10_min_true_likelihood(50), -4.0);
    assert_eq!(java_log10_min_true_likelihood(51), -8.0);
    assert_eq!(java_log10_min_true_likelihood(100), -8.0);
    assert_eq!(java_log10_min_true_likelihood(250), -8.0);
}

#[test]
fn forensic_6r83_java_keeps_equality_drops_strict_less() {
    let thresh = java_log10_min_true_likelihood(100);
    assert_eq!(thresh, -8.0);
    assert!(java_keep_read(-8.0, 100));
    assert!(!java_keep_read(-8.0 - 1e-15, 100));
    assert!(java_keep_read(-7.999, 100));
}

#[test]
fn forensic_6r83_1e6_kernel_noise_cannot_cross_minus_eight() {
    // 6R.82 aligned max |delta| ≈ 7.41e-6. A read whose max sits 1e-4 from the
    // static boundary cannot flip keep/drop from kernel noise.
    let java = -8.0001;
    let rust = java + 7.41e-6;
    assert_eq!(java_keep_read(java, 100), java_keep_read(rust, 100));
}

#[test]
fn forensic_6r83_dump_replay_is_fingerprint_not_raw_index() {
    let a = ReadFp {
        bases: "ACGT".into(),
        bq: "IIII".into(),
        iq: "IIII".into(),
        dq: "IIII".into(),
        gcp: "++++".into(),
    };
    let b = ReadFp {
        bases: "TGCA".into(),
        bq: "IIII".into(),
        iq: "IIII".into(),
        dq: "IIII".into(),
        gcp: "++++".into(),
    };
    let java_order = [&a, &b];
    let rust_order = [&b, &a];
    assert_ne!(java_order, rust_order);
    let js: HashSet<_> = java_order.into_iter().collect();
    let rs: HashSet<_> = rust_order.into_iter().collect();
    assert_eq!(js, rs);
}

#[test]
fn live_post_kernel_read_population_java_vs_rust() {
    if std::env::var("HOLDOUT_6R83").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R83=1");
        return;
    }
    use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
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

    let rust_dump_path = std::env::temp_dir().join("6r83_rust_pairhmm_inputs.txt");
    std::env::set_var(
        "GATK_RS_PAIRHMM_INPUT_DUMP",
        rust_dump_path.to_string_lossy().as_ref(),
    );
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
    let rust_text = std::fs::read_to_string(&rust_dump_path).expect("rust dump");
    let java_regions = parse_dump(&java_text);
    let rust_regions = parse_dump(&rust_text);
    let java_reg = region_with_motif(&java_regions, MOTIF).expect("java motif");
    let rust_reg = region_with_motif(&rust_regions, MOTIF).expect("rust motif");

    let java_kept_i = java_kept_indices(java_reg);
    let java_kept: HashSet<&ReadFp> = java_kept_i.iter().map(|&i| &java_reg.reads[i]).collect();

    let prod_idx: HashSet<usize> = outcome
        .read_likelihoods
        .iter()
        .map(|e| e.read_index.get())
        .collect();
    let prod_haps: HashSet<usize> = outcome
        .read_likelihoods
        .iter()
        .map(|e| e.haplotype_index.get())
        .collect();
    let mut prod_kept: HashSet<&ReadFp> = HashSet::new();
    let mut index_mismatch = 0usize;
    for &i in &prod_idx {
        if i < rust_reg.reads.len() {
            prod_kept.insert(&rust_reg.reads[i]);
        } else {
            index_mismatch += 1;
        }
    }

    let common_prod = java_kept.intersection(&prod_kept).count();
    let java_only_prod: Vec<_> = java_kept.difference(&prod_kept).copied().collect();
    let rust_only_prod: Vec<_> = prod_kept.difference(&java_kept).copied().collect();

    eprintln!(
        "6R.83 PAIRHMM java={}x{} rust={}x{}",
        java_reg.reads.len(),
        java_reg.haps.len(),
        rust_reg.reads.len(),
        rust_reg.haps.len()
    );
    eprintln!(
        "6R.83 java_poorly_modeled_keep={} (Rust kernel dump has no likelihood column; replay uses Java dump)",
        java_kept.len()
    );
    eprintln!(
        "6R.83 production ll_reads={} ll_haps={} geno_reads={} assembly_haps={} index_past_kernel={}",
        prod_idx.len(),
        prod_haps.len(),
        outcome.genotyping_reads.len(),
        outcome.assembly.haplotypes.len(),
        index_mismatch
    );
    eprintln!(
        "6R.83 prod_vs_java_replay COMMON={} JAVA_ONLY={} RUST_ONLY={}",
        common_prod,
        java_only_prod.len(),
        rust_only_prod.len()
    );

    let java_maxs = per_read_max(java_reg);
    let rust_maxs = per_read_max(rust_reg);
    let rust_row_by_fp: HashMap<&ReadFp, usize> = rust_reg
        .reads
        .iter()
        .enumerate()
        .map(|(i, r)| (r, i))
        .collect();

    let mut first_java_only = java_only_prod;
    first_java_only.sort_by_key(|fp| {
        java_reg
            .reads
            .iter()
            .position(|r| r == *fp)
            .unwrap_or(usize::MAX)
    });
    let mut first_rust_only = rust_only_prod;
    first_rust_only.sort_by_key(|fp| {
        rust_reg
            .reads
            .iter()
            .position(|r| r == *fp)
            .unwrap_or(usize::MAX)
    });

    let first_diff = first_java_only
        .first()
        .copied()
        .or_else(|| first_rust_only.first().copied());
    if let Some(fp) = first_diff {
        let ji = java_reg.reads.iter().position(|r| r == fp).unwrap();
        let ri = *rust_row_by_fp.get(fp).unwrap();
        let jmax = java_maxs[ji];
        let rmax = rust_maxs[ri];
        let jthr = java_log10_min_true_likelihood(fp.bq.len().max(1));
        let rec = outcome.genotyping_reads.get(ri).map(|r| r.as_ref());
        let prod_best = outcome
            .read_likelihoods
            .iter()
            .filter(|e| e.read_index.get() == ri)
            .map(|e| e.log10_likelihood)
            .fold(f64::NEG_INFINITY, f64::max);
        let java_ref_lk = java_reg.rows[ji * java_reg.haps.len()].lk;
        let rust_ref_lk = rust_reg.rows[ri * rust_reg.haps.len()].lk;
        eprintln!(
            "6R.83 first_diff_fp bases_hash={:016x} len={} java_row={} rust_row={} java_max={:.10} rust_kernel_max={:.10} prod_best={:.10} thresh={} java_keep={} rust_kernel_keep={} rust_prod_kept={} mapq={:?} start={:?} cigar={:?} unclipped={:?}",
            fnv1a64(&fp.bases),
            fp.bases.len(),
            ji,
            ri,
            jmax,
            rmax,
            prod_best,
            jthr,
            java_keep_read(jmax, fp.bq.len().max(1)),
            java_keep_read(rmax, fp.bq.len().max(1)),
            prod_idx.contains(&ri),
            rec.map(|r| r.mapq()),
            rec.map(|r| r.pos() + 1),
            rec.map(|r| format!("{:?}", r.cigar())),
            rec.map(|r| r.seq().len()),
        );
        eprintln!(
            "6R.83 first_diff hap0_lk java={:?} rust={:?} (column 0; 6R.81 hap_order_equal)",
            java_ref_lk, rust_ref_lk
        );
        let mut near = 0usize;
        for (i, fp) in java_reg.reads.iter().enumerate() {
            let t = java_log10_min_true_likelihood(fp.bq.len().max(1));
            if (java_maxs[i] - t).abs() < 1e-3 {
                near += 1;
            }
        }
        eprintln!("6R.83 reads_within_1e3_of_threshold={near}");
    } else {
        eprintln!("6R.83 first_diff=None (same genotyping fingerprints as Java dump replay)");
    }

    for (label, fps) in [
        ("JAVA_ONLY_PROD", first_java_only.as_slice()),
        ("RUST_ONLY_PROD", first_rust_only.as_slice()),
    ] {
        eprintln!("6R.83 {label} n={}", fps.len());
        for fp in fps.iter().take(20) {
            let ji = java_reg.reads.iter().position(|r| r == *fp);
            let ri = rust_row_by_fp.get(fp).copied();
            let jmax = ji.map(|i| java_maxs[i]);
            let rmax = ri.map(|i| rust_maxs[i]);
            let len = fp.bq.len().max(1);
            eprintln!(
                "  hash={:016x} len={} java_row={:?} rust_row={:?} java_max={:?} rust_max={:?} thresh={} java_keep={:?} rust_kernel_keep={:?} rust_prod={}",
                fnv1a64(&fp.bases),
                fp.bases.len(),
                ji,
                ri,
                jmax,
                rmax,
                java_log10_min_true_likelihood(len),
                jmax.map(|m| java_keep_read(m, len)),
                rmax.map(|m| java_keep_read(m, len)),
                ri.map(|i| prod_idx.contains(&i)).unwrap_or(false),
            );
        }
    }

    let kernel_haps: HashSet<&str> = rust_reg.haps.iter().map(|s| s.as_str()).collect();
    let prod_hap_seqs: HashSet<String> = outcome
        .assembly
        .haplotypes
        .iter()
        .map(|h| String::from_utf8_lossy(&h.bases).into_owned())
        .collect();
    let prod_hap_set: HashSet<&str> = prod_hap_seqs.iter().map(|s| s.as_str()).collect();
    let hap_java_only: Vec<&str> = kernel_haps.difference(&prod_hap_set).copied().collect();
    let hap_rust_only: Vec<&str> = prod_hap_set.difference(&kernel_haps).copied().collect();
    eprintln!(
        "6R.83 haps kernel={} production={} JAVA_ONLY={} RUST_ONLY={}",
        kernel_haps.len(),
        prod_hap_set.len(),
        hap_java_only.len(),
        hap_rust_only.len()
    );
    let ref_pad = outcome
        .assembly
        .haplotypes
        .iter()
        .find(|h| h.is_reference)
        .and_then(|h| h.genome_loc.map(|g| g.start_1based()))
        .unwrap_or(0);
    for s in &hap_java_only {
        let kernel_j = rust_reg.haps.iter().position(|h| h == s);
        let off = POS_SNP.saturating_sub(ref_pad) as usize;
        let b = s.as_bytes().get(off).copied().unwrap_or(b'?');
        eprintln!(
            "6R.83 JAVA_ONLY_HAP hash={:016x} len={} kernel_col={:?} base_at_offset={} ({})",
            fnv1a64(s),
            s.len(),
            kernel_j,
            off,
            b as char
        );
    }
    for s in &hap_rust_only {
        eprintln!(
            "6R.83 RUST_ONLY_HAP hash={:016x} len={}",
            fnv1a64(s),
            s.len()
        );
    }
    let mut n_t = 0usize;
    let mut n_c = 0usize;
    let mut n_other = 0usize;
    for h in &outcome.assembly.haplotypes {
        let pad = h.genome_loc.map(|g| g.start_1based()).unwrap_or(0);
        let off = POS_SNP.saturating_sub(pad) as usize;
        match h.bases.get(off).copied() {
            Some(b'T') => n_t += 1,
            Some(b'C') => n_c += 1,
            _ => n_other += 1,
        }
    }
    eprintln!("6R.83 prod_hap_base_at_offset T={n_t} C={n_c} other={n_other}");

    let live = take_colocated_merge_numerics();
    if let Some(n) = live.iter().find(|n| n.loc == POS_SNP) {
        eprintln!(
            "6R.83 mapper loc={} n_pairhmm_reads={} n_overlap={} n_after_qname={} n_haps={} pools={:?} java_style_pools={:?} merged_ad={:?} subset_ad_remarg={:?} hap_sigs={:?}",
            n.loc,
            n.n_pairhmm_reads,
            n.n_overlap_before_qname_dedupe,
            n.n_reads,
            n.n_haps,
            n.pool_sizes,
            n.java_style_pool_sizes,
            n.merged_ad,
            n.subset_ad_remarginalized,
            n.hap_event_signatures_at_loc
        );
    } else {
        eprintln!("6R.83 mapper numerics missing for canonical loc");
    }

    if let Some(v) = vcf {
        eprintln!(
            "6R.83 rust_vcf pos={} {}:{} PL={:?} AD={:?} QUAL={:?}",
            v.position,
            v.reference,
            v.alternate.join(","),
            v.samples.first().map(|s| &s.pl),
            v.samples.first().map(|s| &s.ad),
            v.quality
        );
    }

    assert_eq!(java_reg.reads.len(), 153);
    assert_eq!(rust_reg.reads.len(), 153);
    assert_eq!(java_reg.haps.len(), 70);
    assert_eq!(rust_reg.haps.len(), 70);
    assert_eq!(java_kept.len(), 136);
    assert_eq!(prod_idx.len(), 136);
    assert_eq!(common_prod, 136);
    assert_eq!(first_java_only.len(), 0);
    assert_eq!(first_rust_only.len(), 0);
    assert_eq!(hap_java_only.len(), 2);
    assert_eq!(hap_rust_only.len(), 0);
    assert_eq!(
        index_mismatch, 0,
        "production LL indices must address kernel rows"
    );
}
