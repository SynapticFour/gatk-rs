//! 6R.82 coordinate-free: PairHMM kernel *matrix* after the common 153×70 population.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`) writes
//! `--pair-hmm-results-file` from `PairHMM.writeToResultsFileIfApplicable`
//! immediately after `subComputeReadLikelihoodGivenHaplotypeLog10`.
//! Layout: hap-bases read-bases BQ IQ DQ GCP expected-result.
//!
//! Compare by **fingerprint**, not raw index:
//!   read  = (bases, BQ, IQ, DQ, GCP)
//!   hap   = sequence
//! Raw-index mismatch with aligned identity is bookkeeping, not a PL cause.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r82_pairhmm_matrix_contract
//! HOLDOUT_6R82=1 cargo test -p gatk-haplotypecaller --test forensic_6r82_pairhmm_matrix_contract live_ -- --nocapture
//! ```

use gatk_haplotypecaller::pairhmm_simd::{
    resolve_pair_hmm_impl, score_read_haps_logless, PairHmmBackend, PairHmmImpl,
};
use gatk_haplotypecaller::{log10_pairhmm_likelihood_exact, logless_pairhmm_likelihood};
use std::collections::{HashMap, HashSet};

#[derive(Clone, Debug, Eq, PartialEq, Hash)]
struct ReadFp {
    bases: String,
    bq: String,
    iq: String,
    dq: String,
    gcp: String,
}

impl ReadFp {
    fn len(&self) -> usize {
        self.bases.len()
    }
}

#[derive(Clone, Debug)]
struct DumpRow {
    hap: String,
    read: ReadFp,
    lk: Option<f64>,
}

#[derive(Clone, Debug)]
struct DumpRegion {
    haps: Vec<String>,
    reads: Vec<ReadFp>,
    /// Row-major: `reads.len() * haps.len()`, hap inner.
    rows: Vec<DumpRow>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum InputField {
    HapBases,
    ReadBases,
    Bq,
    Iq,
    Dq,
    Gcp,
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

fn from_fastq(s: &str) -> Vec<u8> {
    s.bytes().map(|b| b.saturating_sub(33)).collect()
}

fn fnv1a64(s: &str) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for b in s.as_bytes() {
        h ^= u64::from(*b);
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

fn score_kernel(backend: PairHmmBackend, read: &ReadFp, haps: &[String]) -> Vec<f64> {
    let bases = read.bases.as_bytes();
    let bq = from_fastq(&read.bq);
    let iq = from_fastq(&read.iq);
    let dq = from_fastq(&read.dq);
    let gcp = from_fastq(&read.gcp);
    let hap_refs: Vec<&[u8]> = haps.iter().map(|h| h.as_bytes()).collect();
    score_read_haps_logless(backend, bases, &bq, &hap_refs, &iq, &dq, &gcp).expect("kernel")
}

fn first_input_diff(a: &DumpRow, b: &DumpRow) -> Option<(InputField, usize, u8, u8)> {
    for (field, ja, ra) in [
        (InputField::HapBases, a.hap.as_bytes(), b.hap.as_bytes()),
        (
            InputField::ReadBases,
            a.read.bases.as_bytes(),
            b.read.bases.as_bytes(),
        ),
        (InputField::Bq, a.read.bq.as_bytes(), b.read.bq.as_bytes()),
        (InputField::Iq, a.read.iq.as_bytes(), b.read.iq.as_bytes()),
        (InputField::Dq, a.read.dq.as_bytes(), b.read.dq.as_bytes()),
        (
            InputField::Gcp,
            a.read.gcp.as_bytes(),
            b.read.gcp.as_bytes(),
        ),
    ] {
        if ja.len() != ra.len() {
            return Some((
                field,
                ja.len().min(ra.len()),
                ja.len() as u8,
                ra.len() as u8,
            ));
        }
        if let Some(i) = ja.iter().zip(ra.iter()).position(|(x, y)| x != y) {
            return Some((field, i, ja[i], ra[i]));
        }
    }
    None
}

struct AlignedCmp {
    n_cells: usize,
    raw_diff: usize,
    aligned_diff: usize,
    aligned_max_abs: f64,
    aligned_sum_abs: f64,
    first: Option<FirstCell>,
    first_input: Option<(usize, usize, InputField, usize, u8, u8)>,
}

#[derive(Clone, Debug)]
struct FirstCell {
    java_row: usize,
    rust_row: usize,
    hap_j: usize,
    hap_hash: u64,
    hap_len: usize,
    java_lk: f64,
    rust_lk: f64,
    read_len: usize,
    read_hash: u64,
}

fn compare_regions(java: &DumpRegion, rust: &DumpRegion, rust_lks: &[f64]) -> AlignedCmp {
    let jn = java.reads.len() * java.haps.len();
    let rn = rust.reads.len() * rust.haps.len();
    assert_eq!(rust_lks.len(), rn);
    let mut raw_diff = 0usize;
    let n_raw = jn.min(rn);
    for i in 0..n_raw {
        let jl = java.rows[i].lk.unwrap_or(f64::NAN);
        if jl != rust_lks[i] {
            raw_diff += 1;
        }
    }
    raw_diff += jn.abs_diff(rn);

    let rust_row_by_read: HashMap<&ReadFp, usize> =
        rust.reads.iter().enumerate().map(|(i, r)| (r, i)).collect();
    let rust_col_by_hap: HashMap<&str, usize> = rust
        .haps
        .iter()
        .enumerate()
        .map(|(i, h)| (h.as_str(), i))
        .collect();

    let mut aligned_diff = 0usize;
    let mut aligned_max_abs = 0.0f64;
    let mut aligned_sum_abs = 0.0f64;
    let mut first = None;
    let mut first_input = None;
    let mut n_cells = 0usize;
    for (jr, jread) in java.reads.iter().enumerate() {
        let Some(&rr) = rust_row_by_read.get(jread) else {
            continue;
        };
        for (jh, jhap) in java.haps.iter().enumerate() {
            let Some(&rh) = rust_col_by_hap.get(jhap.as_str()) else {
                continue;
            };
            n_cells += 1;
            let ji = jr * java.haps.len() + jh;
            let ri = rr * rust.haps.len() + rh;
            if first_input.is_none() {
                if let Some((field, idx, ja, ra)) = first_input_diff(&java.rows[ji], &rust.rows[ri])
                {
                    first_input = Some((jr, jh, field, idx, ja, ra));
                }
            }
            let jl = java.rows[ji].lk.unwrap_or(f64::NAN);
            let rl = rust_lks[ri];
            let d = (jl - rl).abs();
            aligned_sum_abs += d;
            if d > aligned_max_abs {
                aligned_max_abs = d;
            }
            if jl != rl {
                aligned_diff += 1;
                if first.is_none() {
                    first = Some(FirstCell {
                        java_row: jr,
                        rust_row: rr,
                        hap_j: jh,
                        hap_hash: fnv1a64(jhap),
                        hap_len: jhap.len(),
                        java_lk: jl,
                        rust_lk: rl,
                        read_len: jread.len(),
                        read_hash: fnv1a64(&jread.bases),
                    });
                }
            }
        }
    }
    AlignedCmp {
        n_cells,
        raw_diff,
        aligned_diff,
        aligned_max_abs,
        aligned_sum_abs,
        first,
        first_input,
    }
}

/// Raw index comparison is not the PL comparison when rows are a permutation.
#[test]
fn forensic_6r82_permutation_is_non_causal_when_aligned_cells_match() {
    let haps = vec!["AAA".into(), "AAC".into()];
    let r0 = ReadFp {
        bases: "AA".into(),
        bq: "II".into(),
        iq: "II".into(),
        dq: "II".into(),
        gcp: "++".into(),
    };
    let r1 = ReadFp {
        bases: "AC".into(),
        bq: "II".into(),
        iq: "II".into(),
        dq: "II".into(),
        gcp: "++".into(),
    };
    fn region(haps: &[String], reads: &[ReadFp], lks: &[f64]) -> DumpRegion {
        let mut rows = Vec::new();
        for (ri, read) in reads.iter().enumerate() {
            for (hi, hap) in haps.iter().enumerate() {
                rows.push(DumpRow {
                    hap: hap.clone(),
                    read: read.clone(),
                    lk: Some(lks[ri * haps.len() + hi]),
                });
            }
        }
        DumpRegion {
            haps: haps.to_vec(),
            reads: reads.to_vec(),
            rows,
        }
    }
    let java = region(&haps, &[r0.clone(), r1.clone()], &[-1.0, -2.0, -3.0, -4.0]);
    let rust = region(&haps, &[r1.clone(), r0.clone()], &[-3.0, -4.0, -1.0, -2.0]);
    let rust_lks = vec![-3.0, -4.0, -1.0, -2.0];
    let cmp = compare_regions(&java, &rust, &rust_lks);
    assert!(
        cmp.raw_diff > 0,
        "permuted rows must fail raw index compare"
    );
    assert_eq!(cmp.aligned_diff, 0);
    assert_eq!(cmp.aligned_max_abs, 0.0);
    assert!(cmp.first.is_none());
    assert!(cmp.first_input.is_none());
}

/// Java-order row-major scan reports the first mismatched aligned cell.
#[test]
fn forensic_6r82_first_divergent_cell_is_java_row_major() {
    let haps: Vec<String> = vec!["AAA".into(), "AAC".into()];
    let r0 = ReadFp {
        bases: "AA".into(),
        bq: "II".into(),
        iq: "II".into(),
        dq: "II".into(),
        gcp: "++".into(),
    };
    let r1 = ReadFp {
        bases: "AC".into(),
        bq: "II".into(),
        iq: "II".into(),
        dq: "II".into(),
        gcp: "++".into(),
    };
    let java_rows = vec![
        DumpRow {
            hap: haps[0].clone(),
            read: r0.clone(),
            lk: Some(-1.0),
        },
        DumpRow {
            hap: haps[1].clone(),
            read: r0.clone(),
            lk: Some(-2.0),
        },
        DumpRow {
            hap: haps[0].clone(),
            read: r1.clone(),
            lk: Some(-3.0),
        },
        DumpRow {
            hap: haps[1].clone(),
            read: r1.clone(),
            lk: Some(-4.0),
        },
    ];
    let java = DumpRegion {
        haps: haps.clone(),
        reads: vec![r0.clone(), r1.clone()],
        rows: java_rows,
    };
    let rust = java.clone();
    let rust_lks = vec![-1.0, -2.0, -3.0, -9.0];
    let cmp = compare_regions(&java, &rust, &rust_lks);
    assert_eq!(cmp.aligned_diff, 1);
    let first = cmp.first.expect("cell");
    assert_eq!(first.java_row, 1);
    assert_eq!(first.hap_j, 1);
    assert_eq!(first.java_lk, -4.0);
    assert_eq!(first.rust_lk, -9.0);
}

#[test]
fn forensic_6r82_production_backend_is_resolved_not_guessed() {
    let backend = resolve_pair_hmm_impl(PairHmmImpl::FastestAvailable);
    eprintln!("6R.82 production FastestAvailable -> {}", backend.label());
    assert!(
        matches!(
            backend,
            PairHmmBackend::Avx2F64
                | PairHmmBackend::NeonF64
                | PairHmmBackend::PackedF64
                | PairHmmBackend::LoglessScalar
        ),
        "unexpected FastestAvailable backend {}",
        backend.label()
    );
}

/// Java `--pair-hmm-results-file` prints `String.format("%e", lk)` (~6 decimal digits).
/// Exact f64 inequality against that dump is not a material kernel miss.
#[test]
fn forensic_6r82_java_percent_e_dump_rounding_is_not_a_kernel_miss() {
    // Java `String.format("%e", x)` keeps ~6 digits after the decimal.
    let java_dump = -9.806065f64;
    let rust_full = -9.806065085111925f64;
    assert_ne!(java_dump, rust_full);
    let d = (java_dump - rust_full).abs();
    assert!(
        d < 1e-6,
        "percent-e rounding must stay at dump ULP, not PL scale: {d}"
    );
}

#[test]
fn forensic_6r82_kernel_families_run_on_identical_input_tuple() {
    let read = b"ACGTACGT";
    let hap = b"ACGTACGTAC";
    let bq = vec![30u8; read.len()];
    let gop = vec![45u8; read.len()];
    let gcp = vec![10u8; read.len()];
    let log10 = log10_pairhmm_likelihood_exact(read, &bq, hap, &gop, &gop, &gcp).unwrap();
    let logless = logless_pairhmm_likelihood(read, &bq, hap, &gop, &gop, &gcp).unwrap();
    assert!(log10 <= 0.0 && log10.is_finite());
    assert!(logless <= 0.0 && logless.is_finite());
    eprintln!(
        "6R.82 scalar log10={log10} logless={logless} abs_delta={}",
        (log10 - logless).abs()
    );
}

#[test]
fn live_pairhmm_matrix_java_vs_rust() {
    if std::env::var("HOLDOUT_6R82").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R82=1");
        return;
    }
    use gatk_core::reference::{parse_intervals_cli_string, SequenceDictionary};
    use gatk_haplotypecaller::{
        call_disposition, flatten_assembly_regions, traverse_assembly_region_walker,
        try_emit_call_region_variants, AssemblyRegionCallDisposition, CallRegionArgs,
        HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
        DEFAULT_STAND_EMIT_CONFIDENCE,
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

    let rust_dump_path = std::env::temp_dir().join("6r82_rust_pairhmm_inputs.txt");
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
    let java_reg = region_with_motif(&java_regions, MOTIF).expect("java motif region");
    let rust_reg = region_with_motif(&rust_regions, MOTIF).expect("rust motif region");

    let java_reads: HashSet<&ReadFp> = java_reg.reads.iter().collect();
    let rust_reads: HashSet<&ReadFp> = rust_reg.reads.iter().collect();
    let java_only_reads = java_reads.difference(&rust_reads).count();
    let rust_only_reads = rust_reads.difference(&java_reads).count();
    let common_reads = java_reads.intersection(&rust_reads).count();
    let hap_order_equal = java_reg.haps == rust_reg.haps;
    let read_order_equal = java_reg.reads == rust_reg.reads;
    let java_dup = {
        let mut c = HashMap::<&ReadFp, usize>::new();
        for r in &java_reg.reads {
            *c.entry(r).or_default() += 1;
        }
        c.values().filter(|n| **n > 1).count()
    };
    let rust_dup = {
        let mut c = HashMap::<&ReadFp, usize>::new();
        for r in &rust_reg.reads {
            *c.entry(r).or_default() += 1;
        }
        c.values().filter(|n| **n > 1).count()
    };

    let backend = resolve_pair_hmm_impl(PairHmmImpl::FastestAvailable);
    let mut rust_prod = Vec::with_capacity(rust_reg.reads.len() * rust_reg.haps.len());
    let mut rust_log10 = Vec::with_capacity(rust_prod.capacity());
    let mut rust_logless = Vec::with_capacity(rust_prod.capacity());
    for read in &rust_reg.reads {
        rust_prod.extend(score_kernel(backend, read, &rust_reg.haps));
        rust_log10.extend(score_kernel(
            PairHmmBackend::Log10Scalar,
            read,
            &rust_reg.haps,
        ));
        rust_logless.extend(score_kernel(
            PairHmmBackend::LoglessScalar,
            read,
            &rust_reg.haps,
        ));
    }
    let mut java_inputs_prod = Vec::with_capacity(java_reg.reads.len() * java_reg.haps.len());
    for read in &java_reg.reads {
        java_inputs_prod.extend(score_kernel(backend, read, &java_reg.haps));
    }

    let cmp_prod = compare_regions(java_reg, rust_reg, &rust_prod);
    let cmp_log10 = compare_regions(java_reg, rust_reg, &rust_log10);
    let cmp_logless = compare_regions(java_reg, rust_reg, &rust_logless);
    let mut kernel_on_java_diff = 0usize;
    let mut kernel_on_java_max = 0.0f64;
    for (i, row) in java_reg.rows.iter().enumerate() {
        let jl = row.lk.unwrap_or(f64::NAN);
        let rl = java_inputs_prod[i];
        let d = (jl - rl).abs();
        if jl != rl {
            kernel_on_java_diff += 1;
        }
        if d > kernel_on_java_max {
            kernel_on_java_max = d;
        }
    }

    eprintln!(
        "6R.82 dims java={}x{} rust={}x{} backend={}",
        java_reg.reads.len(),
        java_reg.haps.len(),
        rust_reg.reads.len(),
        rust_reg.haps.len(),
        backend.label()
    );
    eprintln!(
        "6R.82 hap_order_equal={hap_order_equal} hap_set_equal={}",
        java_reg.haps.iter().collect::<HashSet<_>>()
            == rust_reg.haps.iter().collect::<HashSet<_>>()
    );
    eprintln!(
        "6R.82 reads COMMON={common_reads} JAVA_ONLY={java_only_reads} RUST_ONLY={rust_only_reads} dup_java={java_dup} dup_rust={rust_dup} read_order_equal={read_order_equal}"
    );
    eprintln!(
        "6R.82 raw_diff={} aligned_prod diff={} max_abs={:.6e} mean_abs={:.6e} n_cells={}",
        cmp_prod.raw_diff,
        cmp_prod.aligned_diff,
        cmp_prod.aligned_max_abs,
        if cmp_prod.n_cells == 0 {
            0.0
        } else {
            cmp_prod.aligned_sum_abs / cmp_prod.n_cells as f64
        },
        cmp_prod.n_cells
    );
    eprintln!(
        "6R.82 aligned_log10_diff={} max={:.6e} aligned_logless_diff={} max={:.6e}",
        cmp_log10.aligned_diff,
        cmp_log10.aligned_max_abs,
        cmp_logless.aligned_diff,
        cmp_logless.aligned_max_abs
    );
    eprintln!(
        "6R.82 rust_kernel_on_java_inputs diff={} max_abs={:.6e} (isolates kernel given Java tuples)",
        kernel_on_java_diff, kernel_on_java_max
    );
    eprintln!(
        "6R.82 first_input={:?} (None = identical kernel tuples)",
        cmp_prod.first_input
    );
    eprintln!("6R.82 first_cell={:?}", cmp_prod.first);
    if let Some(c) = &cmp_prod.first {
        eprintln!(
            "6R.82 first_cell java_row={} rust_row={} hap_j={} hap_hash={:016x} hap_len={} read_len={} read_bases_hash={:016x} java_lk={} rust_lk={} delta={:.6e}",
            c.java_row,
            c.rust_row,
            c.hap_j,
            c.hap_hash,
            c.hap_len,
            c.read_len,
            c.read_hash,
            c.java_lk,
            c.rust_lk,
            (c.java_lk - c.rust_lk).abs()
        );
    }
    if let Some(v) = vcf {
        eprintln!(
            "6R.82 rust_vcf pos={} {}:{} PL={:?} AD={:?} QUAL={:?}",
            v.position,
            v.reference,
            v.alternate.join(","),
            v.samples.first().map(|s| &s.pl),
            v.samples.first().map(|s| &s.ad),
            v.quality
        );
    }

    assert_eq!(java_reg.reads.len(), 153);
    assert_eq!(java_reg.haps.len(), 70);
    assert_eq!(rust_reg.reads.len(), 153);
    assert_eq!(rust_reg.haps.len(), 70);
    assert!(hap_order_equal, "6R.81: haplotype order must remain equal");
    assert_eq!(java_only_reads, 0);
    assert_eq!(rust_only_reads, 0);
    assert_eq!(common_reads, 153);
    assert_eq!(java_dup, 0);
    assert_eq!(rust_dup, 0);
    assert_eq!(cmp_prod.n_cells, 153 * 70);
    assert!(
        cmp_prod.first_input.is_none(),
        "aligned PairHMM input tuples (bases/BQ/IQ/DQ/GCP/hap) must match; got {:?}",
        cmp_prod.first_input
    );
    assert!(
        cmp_prod.aligned_max_abs < 1e-5,
        "aligned kernel |delta| must stay at Java %e dump / f32 noise scale, max={}",
        cmp_prod.aligned_max_abs
    );
    assert!(
        kernel_on_java_max < 1e-5,
        "Rust kernel on Java input tuples must match dump at %e scale, max={kernel_on_java_max}"
    );
}
