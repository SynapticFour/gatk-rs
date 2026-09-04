//! 6R.76 coordinate-free: PairHMM haplotype sequence SET.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//! `AssemblyResultSet.trimDownHaplotypes` uses `h.trim(span, true)` then
//! HashMap equality on post-trim bases (uniqueness defaults to 0).
//! PairHMM haplotypes are `assemblyResult.getHaplotypeList()` **before**
//! `filterAlleles`.
//!
//! Live kernel dump (SNP-motif region) after 6R.75:
//! `only_java = 2` (J0, J1), `only_rust = 1` (R0). J0 vs R0 is Hamming-1 at
//! index 226 (`C` vs `G`). That difference already exists on the untrimmed
//! k=25 `407M` path — trim/padding/dedup cannot invent J0 from R0.
//!
//! No production change in this round (first divergence is assembly/graph
//! sequence, frozen). Case E (`ignoreRefState=false` keeping a ref+alt
//! duplicate) is the 6R.81 PairHMM-population arrow, not stacked here.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r76_haplotype_set_contract
//! HOLDOUT_6R76=1 cargo test -p gatk-haplotypecaller --test forensic_6r76_haplotype_set_contract live_haplotype -- --nocapture
//! ```

use gatk_haplotypecaller::haplotype::haplotype_size_and_base_order;
use gatk_haplotypecaller::{
    AssemblyRegion, AssemblyResultSet, AssemblyStatus, Cigar, CigarOperator, FeatureContext,
    GenomeLoc, GenomePosition, Haplotype, ReferenceContext,
};
use std::collections::{HashMap, HashSet};

/// Java-only PairHMM sequences from the 6R.75 kernel dump (SNP-motif region).
const JAVA_ONLY_J0: &[u8] = b"CATGGAGCCTGACCTTATTTGAAGTAGGGCATTTGCAGATGTATTTAAGATATTTGAGGCTGGGCACAGTGGCTCACGTCTGTAATCCCAGCACTTTGAAAGGCCGAGGCAGGTGGATTCACCTGAGGTCAGGAGTTTGAGACCAGCCTGTCCCACATGGTGAAAAGCCCGTATCTACCAAAAATACAAACGTTAGCTGTGTGTGGTGGTGGCGGCACCTGTAATCCCAGCTACTCGAGAGCCAGAG";
const JAVA_ONLY_J1: &[u8] = b"CATGGAGCCTGACTTTATTTGAAGTAGGGCATTTGCAGATGTATTTAAGATATTTGAGGCCGGGCACAGTGGCTCATGTCTGTAATCCCAGCACTTTAAAAGGCTGAGGCAGGTGTATTCACCTGAGGTCAGGAGTTCGAGACCAGCCTGGCCAACATGGTGAAAAGCCCGTATCTACCAAAAATACAAAAGTTAGCTGGGTGTGGTGGTGGCGGCACCTGTAATCCCAGCTACTCGAGAGCCAGAG";
const RUST_ONLY_R0: &[u8] = b"CATGGAGCCTGACCTTATTTGAAGTAGGGCATTTGCAGATGTATTTAAGATATTTGAGGCTGGGCACAGTGGCTCACGTCTGTAATCCCAGCACTTTGAAAGGCCGAGGCAGGTGGATTCACCTGAGGTCAGGAGTTTGAGACCAGCCTGTCCCACATGGTGAAAAGCCCGTATCTACCAAAAATACAAACGTTAGCTGTGTGTGGTGGTGGCGGCACCTGTAATCGCAGCTACTCGAGAGCCAGAG";

const UNTRIMMED_PREFIX: usize = 74;
const UNTRIMMED_SUFFIX: usize = 86;

/// GATK 4.4.0.0 `AssemblyResultSet.trimDownHaplotypes` replay:
/// `h.trim(span, true)` then HashMap keyed by post-trim bases (all trimmed
/// haplotypes are non-ref, uniqueness 0), restoring the reference flag from
/// the original haplotype that won the collision.
fn java_trim_down_haplotypes(haps: &[Haplotype], span: &GenomeLoc) -> Vec<Haplotype> {
    let mut original_by_trimmed: Vec<(Haplotype, usize)> = Vec::new();
    let mut index_by_bases: HashMap<Vec<u8>, usize> = HashMap::new();
    for (orig_i, h) in haps.iter().enumerate() {
        let Some(trimmed) = h.trim(span, true) else {
            continue;
        };
        if let Some(&idx) = index_by_bases.get(&trimmed.bases) {
            if h.is_reference {
                original_by_trimmed[idx] = (trimmed, orig_i);
            }
        } else {
            index_by_bases.insert(trimmed.bases.clone(), original_by_trimmed.len());
            original_by_trimmed.push((trimmed, orig_i));
        }
    }
    let mut out = Vec::with_capacity(original_by_trimmed.len());
    for (mut trimmed, orig_i) in original_by_trimmed {
        if haps[orig_i].is_reference {
            trimmed.is_reference = true;
        }
        out.push(trimmed);
    }
    out.sort_by(haplotype_size_and_base_order);
    out
}

fn hap_match(bases: &[u8], is_ref: bool, loc: GenomeLoc) -> Haplotype {
    let mut h = Haplotype::new(bases.to_vec(), is_ref);
    let mut c = Cigar::new();
    c.push(bases.len(), CigarOperator::Match);
    h.cigar = Some(c);
    h.genome_loc = Some(loc);
    h
}

fn dummy_region(start: u64, end: u64) -> AssemblyRegion {
    AssemblyRegion {
        contig: "synth".into(),
        start: GenomePosition::new_1based(start),
        end: GenomePosition::new_1based(end),
        is_active: true,
        extended_start: GenomePosition::new_1based(start),
        extended_end: GenomePosition::new_1based(end),
        extension: 0,
        reads: Vec::new(),
        read_qnames: Vec::new(),
        reference: ReferenceContext::empty(),
        features: FeatureContext::empty(),
        pileup_loci: Vec::new(),
    }
}

fn r0_untrimmed_bases() -> Vec<u8> {
    let mut v = vec![b'A'; UNTRIMMED_PREFIX];
    v.extend_from_slice(RUST_ONLY_R0);
    v.extend_from_slice(&vec![b'T'; UNTRIMMED_SUFFIX]);
    v
}

fn unique_seqs(haps: &[Haplotype]) -> HashSet<Vec<u8>> {
    haps.iter().map(|h| h.bases.clone()).collect()
}

fn hamming(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b.iter()).filter(|(x, y)| x != y).count() + a.len().abs_diff(b.len())
}

fn label_known(seq: &[u8]) -> &'static str {
    if seq == JAVA_ONLY_J0 {
        "J0"
    } else if seq == JAVA_ONLY_J1 {
        "J1"
    } else if seq == RUST_ONLY_R0 {
        "R0"
    } else {
        "-"
    }
}

fn substring_hits(untrimmed: &[Haplotype], needle: &[u8]) -> Vec<usize> {
    untrimmed
        .iter()
        .enumerate()
        .filter(|(_, h)| h.bases.windows(needle.len()).any(|w| w == needle))
        .map(|(i, _)| i)
        .collect()
}

fn snp_motif_haps(dump: &str, motif: &str) -> Vec<String> {
    let mut parsed: Vec<(String, String)> = Vec::new();
    for line in dump.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut p = line.split(' ');
        let hap = p.next().unwrap_or("").to_string();
        let read = p.next().unwrap_or("").to_string();
        parsed.push((hap, read));
    }
    let mut i = 0;
    while i < parsed.len() {
        let rid = parsed[i].1.clone();
        let mut haps = Vec::new();
        let mut j = i;
        while j < parsed.len() && parsed[j].1 == rid {
            haps.push(parsed[j].0.clone());
            j += 1;
        }
        if haps.iter().any(|h| h.contains(motif)) {
            return haps;
        }
        i = j;
    }
    Vec::new()
}

#[test]
fn forensic_6r76_java_only_and_rust_only_are_distinct_sequences() {
    assert_ne!(JAVA_ONLY_J0, JAVA_ONLY_J1);
    assert_ne!(JAVA_ONLY_J0, RUST_ONLY_R0);
    assert_ne!(JAVA_ONLY_J1, RUST_ONLY_R0);
    assert_eq!(JAVA_ONLY_J0.len(), 247);
    assert_eq!(JAVA_ONLY_J1.len(), 247);
    assert_eq!(RUST_ONLY_R0.len(), 247);
    let diffs: Vec<usize> = JAVA_ONLY_J0
        .iter()
        .zip(RUST_ONLY_R0.iter())
        .enumerate()
        .filter(|(_, (a, b))| a != b)
        .map(|(i, _)| i)
        .collect();
    assert_eq!(diffs, vec![226], "J0 vs R0 must remain a single SNP at 226");
    assert_eq!(JAVA_ONLY_J0[226], b'C');
    assert_eq!(RUST_ONLY_R0[226], b'G');
    assert_ne!(
        RUST_ONLY_R0, JAVA_ONLY_J0,
        "C is not A after any trim/pad interpretation"
    );
    assert_ne!(
        RUST_ONLY_R0, JAVA_ONLY_J1,
        "C is not B after any trim/pad interpretation"
    );
}

#[test]
fn forensic_6r76_java_trim_down_keeps_j0_and_r0_as_two_sequences() {
    let loc = GenomeLoc::new(1, 247);
    let haps = vec![
        hap_match(JAVA_ONLY_J0, false, loc),
        hap_match(RUST_ONLY_R0, false, loc),
        hap_match(JAVA_ONLY_J1, false, loc),
    ];
    let trimmed = java_trim_down_haplotypes(&haps, &loc);
    let set = unique_seqs(&trimmed);
    assert!(set.contains(JAVA_ONLY_J0), "J0 must survive Java trimDown");
    assert!(set.contains(RUST_ONLY_R0), "R0 must survive Java trimDown");
    assert!(set.contains(JAVA_ONLY_J1), "J1 must survive Java trimDown");
    assert_eq!(set.len(), 3, "Hamming-1 SNPs are not Java trim duplicates");
}

#[test]
fn forensic_6r76_trim_of_r0_untrimmed_window_cannot_invent_j0() {
    let pad = 1000u64;
    let untrimmed_len = (UNTRIMMED_PREFIX + RUST_ONLY_R0.len() + UNTRIMMED_SUFFIX) as u64;
    let full = GenomeLoc::new(pad, pad + untrimmed_len - 1);
    let trim_start = pad + UNTRIMMED_PREFIX as u64;
    let trim_end = trim_start + RUST_ONLY_R0.len() as u64 - 1;
    let span = GenomeLoc::new(trim_start, trim_end);

    let untrimmed_bases = r0_untrimmed_bases();
    assert_eq!(
        &untrimmed_bases[UNTRIMMED_PREFIX..UNTRIMMED_PREFIX + RUST_ONLY_R0.len()],
        RUST_ONLY_R0
    );

    let mut ref_window = RUST_ONLY_R0.to_vec();
    ref_window[13] = b'T';
    let mut ref_untrimmed = vec![b'A'; UNTRIMMED_PREFIX];
    ref_untrimmed.extend_from_slice(&ref_window);
    ref_untrimmed.extend_from_slice(&vec![b'T'; UNTRIMMED_SUFFIX]);

    let haps = vec![
        hap_match(&untrimmed_bases, false, full),
        hap_match(&ref_untrimmed, true, full),
    ];
    let java_replay = java_trim_down_haplotypes(&haps, &span);
    let java_set = unique_seqs(&java_replay);
    assert!(
        java_set.contains(RUST_ONLY_R0),
        "trim window of the R0-bearing path is R0"
    );
    assert!(
        !java_set.contains(JAVA_ONLY_J0),
        "Java trimDown must not flip offset 226 C←G"
    );
    assert!(
        !java_set.contains(JAVA_ONLY_J1),
        "Java trimDown must not invent J1 from an R0 path"
    );

    let assembly = AssemblyResultSet::from_assembly_for_calling_owned(
        AssemblyStatus::AssembledSomeVariation,
        25,
        haps,
        ref_untrimmed.as_slice(),
        pad,
        "synth",
        0,
    );
    let region = dummy_region(trim_start, trim_end);
    let rust_trim = assembly.trim_to(&region).expect("trim_to");
    let rust_set = unique_seqs(&rust_trim.haplotypes);
    assert!(
        rust_set.contains(RUST_ONLY_R0),
        "production trim_to must emit the R0 window"
    );
    assert!(
        !rust_set.contains(JAVA_ONLY_J0),
        "production trim_to must not invent J0 (pre-6R.76 assembly divergence)"
    );
}

/// Java `trimDownHaplotypes` note: ref and non-ref that trim to the same bases
/// collapse, and the survivor is marked reference. 6R.81 production `trim_to`
/// now matches this; kept here as the 6R.76 Java-replay pin.
#[test]
fn forensic_6r76_java_trim_down_collapses_ref_and_alt_identical_bases() {
    let loc = GenomeLoc::new(1, 20);
    let bases = b"ACGTACGTACGTACGTACGT";
    let haps = vec![hap_match(bases, false, loc), hap_match(bases, true, loc)];
    let trimmed = java_trim_down_haplotypes(&haps, &loc);
    assert_eq!(trimmed.len(), 1);
    assert!(trimmed[0].is_reference);
    assert_eq!(trimmed[0].bases.as_slice(), bases);
}

/// Live pipeline dump. Not the sole proof.
#[test]
fn live_haplotype_pipeline_stages() {
    if std::env::var("HOLDOUT_6R76").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R76=1");
        return;
    }
    use gatk_core::reference::{
        parse_intervals_cli_string, ReferenceWindowCache, SequenceDictionary,
    };
    use gatk_haplotypecaller::{
        assemble_reads_with_finalized, call_disposition, flatten_assembly_regions,
        traverse_assembly_region_walker, AssemblyRegionCallDisposition, CallRegionArgs,
        HaplotypeCallerEngine, ReadFilterParams, WalkerTraversalConfig,
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
    if !ref_fasta.is_file() || !bam.is_file() {
        eprintln!("skip: live BAM/ref missing");
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
    let region = covering[0];
    let args = CallRegionArgs::strict_java();
    let mut ref_cache = ReferenceWindowCache::new(ref_fasta.clone(), 4);
    let mut owned = region.clone();
    let assembled =
        assemble_reads_with_finalized(&mut owned, &dict, &mut ref_cache, &args.assemble)
            .expect("assemble");
    let untrimmed = assembled.assembly;

    eprintln!(
        "UNTRIMMED n={} unique={} n_ref={}",
        untrimmed.haplotypes.len(),
        unique_seqs(&untrimmed.haplotypes).len(),
        untrimmed
            .haplotypes
            .iter()
            .filter(|h| h.is_reference)
            .count()
    );
    for (name, needle) in [
        ("J0", JAVA_ONLY_J0),
        ("J1", JAVA_ONLY_J1),
        ("R0", RUST_ONLY_R0),
    ] {
        let hits = substring_hits(&untrimmed.haplotypes, needle);
        eprintln!(
            "  untrimmed contains {name} as substring: {} indices={hits:?}",
            !hits.is_empty()
        );
    }
    assert!(
        substring_hits(&untrimmed.haplotypes, JAVA_ONLY_J0).len() == 1,
        "6R.80: live untrimmed assembly retains J0 (Java K=128 rank 124)"
    );
    assert!(
        substring_hits(&untrimmed.haplotypes, RUST_ONLY_R0).is_empty(),
        "6R.80: R0 is the 129th SeqGraph sink and is outside production K=128"
    );

    let rust_dump_path = std::env::temp_dir().join("6r76_rust_pairhmm_inputs.txt");
    std::env::set_var(
        "GATK_RS_PAIRHMM_INPUT_DUMP",
        rust_dump_path.to_string_lossy().as_ref(),
    );
    let outcome = HaplotypeCallerEngine::call_region(region, &dict, &ref_fasta, &args)
        .expect("call")
        .expect("Some");
    let span = outcome
        .assembly
        .haplotypes
        .iter()
        .find_map(|h| h.genome_loc)
        .expect("trimmed genome_loc");

    let mut trimmed_region = region.clone();
    trimmed_region.extended_start = GenomePosition::new_1based(span.start_1based());
    trimmed_region.extended_end = GenomePosition::new_1based(span.end_1based());
    let rust_trim = untrimmed.trim_to(&trimmed_region).expect("trim_to");
    let java_replay = java_trim_down_haplotypes(&untrimmed.haplotypes, &span);
    let rust_set = unique_seqs(&rust_trim.haplotypes);
    let java_replay_set = unique_seqs(&java_replay);
    assert_eq!(
        rust_set, java_replay_set,
        "Java trimDown replay of rust untrimmed must match rust unique trim_to SET"
    );
    assert!(rust_set.contains(JAVA_ONLY_J0));
    assert!(!rust_set.contains(RUST_ONLY_R0));

    if rust_dump_path.is_file() && java_dump_path.is_file() {
        let rust_dump = std::fs::read_to_string(&rust_dump_path).expect("rust dump");
        let java_dump = std::fs::read_to_string(&java_dump_path).expect("java dump");
        let rust_haps = snp_motif_haps(&rust_dump, MOTIF);
        let java_haps = snp_motif_haps(&java_dump, MOTIF);
        let rs: HashSet<&str> = rust_haps.iter().map(|s| s.as_str()).collect();
        let js: HashSet<&str> = java_haps.iter().map(|s| s.as_str()).collect();
        let only_j: Vec<_> = js.difference(&rs).copied().collect();
        let only_r: Vec<_> = rs.difference(&js).copied().collect();
        eprintln!(
            "PAIRHMM dump only_java={} only_rust={}",
            only_j.len(),
            only_r.len()
        );
        for s in &only_j {
            eprintln!("  JAVA ONLY {} len={}", label_known(s.as_bytes()), s.len());
        }
        for s in &only_r {
            eprintln!("  RUST ONLY {} len={}", label_known(s.as_bytes()), s.len());
        }
        assert!(
            rs.iter().any(|s| s.as_bytes() == JAVA_ONLY_J0),
            "6R.80: rust PairHMM dump contains J0"
        );
        assert!(
            !rs.iter().any(|s| s.as_bytes() == RUST_ONLY_R0),
            "6R.80: rust PairHMM dump no longer contains R0"
        );
        assert_eq!(hamming(JAVA_ONLY_J0, RUST_ONLY_R0), 1);
    }
}
