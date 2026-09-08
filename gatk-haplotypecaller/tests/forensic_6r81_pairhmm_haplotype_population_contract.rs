//! 6R.81 coordinate-free: PairHMM haplotype *population* after trimDown.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//! `AssemblyResultSet.trimDownHaplotypes`:
//!   `h.trim(span, true)` then `HashMap` keyed by `Haplotype.equals`.
//! After `ignoreRefState=true`, uniqueness is 0 and every trimmed haplotype is
//! non-ref, so equality is **bases only**. A reference haplotype and an alt
//! that trim to the same bases collapse; the survivor is marked reference.
//!
//! Pre-6R.81 production `trim_to` used `h.trim(span, false)` and keyed
//! `(bases, is_reference)`, which kept a ref+alt duplicate of a common
//! sequence as two PairHMM columns (kernel 153×71 vs Java 153×70) while the
//! unique-sequence SET already matched.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r81_pairhmm_haplotype_population_contract
//! HOLDOUT_6R81=1 cargo test -p gatk-haplotypecaller --test forensic_6r81_pairhmm_haplotype_population_contract live_ -- --nocapture
//! ```

use gatk_haplotypecaller::haplotype::haplotype_size_and_base_order;
use gatk_haplotypecaller::{
    AssemblyRegion, AssemblyResultSet, AssemblyStatus, Cigar, CigarOperator, FeatureContext,
    GenomeLoc, GenomePosition, Haplotype, ReferenceContext,
};
use std::collections::{HashMap, HashSet};

const BASES_A: &[u8] = b"ACGTACGTACGTACGTACGT";
const BASES_B: &[u8] = b"ACGTACGTACGTACGTACGG";

/// Pre-6R.81 production `trim_to`: `trim(span, false)` + key `(bases, is_reference)`.
fn legacy_trim_keep_ref_state(haps: &[Haplotype], span: &GenomeLoc) -> Vec<Haplotype> {
    let mut trimmed_list: Vec<Haplotype> = Vec::new();
    let mut index_by_hap_key: HashMap<(Vec<u8>, bool), usize> = HashMap::new();
    for h in haps {
        let Some(t) = h.trim(span, false) else {
            continue;
        };
        let key = (t.bases.clone(), t.is_reference);
        if let Some(&idx) = index_by_hap_key.get(&key) {
            if h.is_reference {
                trimmed_list[idx] = t;
            }
        } else {
            index_by_hap_key.insert(key, trimmed_list.len());
            trimmed_list.push(t);
        }
    }
    trimmed_list.sort_by(haplotype_size_and_base_order);
    trimmed_list
}

/// GATK 4.4.0.0 `trimDownHaplotypes` replay (same as 6R.76 helper).
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

fn production_trim_to(
    haps: Vec<Haplotype>,
    ref_bases: &[u8],
    start: u64,
    end: u64,
) -> Vec<Haplotype> {
    let assembly = AssemblyResultSet::from_assembly_for_calling_owned(
        AssemblyStatus::AssembledSomeVariation,
        25,
        haps,
        ref_bases,
        start,
        "synth",
        0,
    );
    let region = dummy_region(start, end);
    assembly.trim_to(&region).expect("trim_to").haplotypes
}

fn unique_seqs(haps: &[Haplotype]) -> HashSet<Vec<u8>> {
    haps.iter().map(|h| h.bases.clone()).collect()
}

fn fnv1a64(seq: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in seq {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

/// Java `trimDownHaplotypes` + production `trim_to` collapse a ref+alt copy of
/// the same bases to one haplotype marked reference. Legacy `(bases, is_ref)`
/// keying keeps both PairHMM columns.
#[test]
fn forensic_6r81_ref_alt_identical_bases_collapse_to_one_reference() {
    let loc = GenomeLoc::new(1, BASES_A.len() as u64);
    let haps = vec![
        hap_match(BASES_A, false, loc),
        hap_match(BASES_A, true, loc),
    ];
    let legacy = legacy_trim_keep_ref_state(&haps, &loc);
    let java = java_trim_down_haplotypes(&haps, &loc);
    let prod = production_trim_to(haps.clone(), BASES_A, 1, BASES_A.len() as u64);

    assert_eq!(
        legacy.len(),
        2,
        "pre-6R.81 policy keeps ref+alt as two objects"
    );
    assert_eq!(unique_seqs(&legacy).len(), 1);
    assert_eq!(
        legacy.iter().filter(|h| h.is_reference).count(),
        1,
        "legacy retains one ref and one alt of the same bases"
    );
    assert_eq!(legacy.iter().filter(|h| !h.is_reference).count(), 1);

    assert_eq!(java.len(), 1);
    assert!(java[0].is_reference);
    assert_eq!(java[0].bases.as_slice(), BASES_A);

    assert_eq!(
        prod.len(),
        1,
        "production trim_to must match Java object count"
    );
    assert!(prod[0].is_reference);
    assert_eq!(prod[0].bases.as_slice(), BASES_A);
    assert_eq!(unique_seqs(&prod).len(), 1);
}

/// Collision order: reference arriving after the alt still marks the survivor reference.
#[test]
fn forensic_6r81_later_reference_replaces_alt_and_keeps_ref_flag() {
    let loc = GenomeLoc::new(1, BASES_A.len() as u64);
    let haps = vec![
        hap_match(BASES_A, true, loc),
        hap_match(BASES_A, false, loc),
    ];
    let prod = production_trim_to(haps, BASES_A, 1, BASES_A.len() as u64);
    assert_eq!(prod.len(), 1);
    assert!(prod[0].is_reference);
}

/// Two alts with identical bases also collapse (Java HashMap after trim(true)).
/// A distinct reference haplotype is required so `tag_padded_reference_span`
/// can assign a valid genome loc (Java `trimTo` refuses a missing ref hap).
#[test]
fn forensic_6r81_two_alts_identical_bases_collapse() {
    let loc = GenomeLoc::new(1, BASES_A.len() as u64);
    let haps = vec![
        hap_match(BASES_B, true, loc),
        hap_match(BASES_A, false, loc),
        hap_match(BASES_A, false, loc),
    ];
    let prod = production_trim_to(haps.clone(), BASES_B, 1, BASES_B.len() as u64);
    let java = java_trim_down_haplotypes(&haps, &loc);
    assert_eq!(java.len(), 2);
    assert_eq!(prod.len(), 2);
    assert_eq!(unique_seqs(&prod).len(), 2);
    assert!(prod
        .iter()
        .any(|h| h.is_reference && h.bases.as_slice() == BASES_B));
    assert_eq!(
        prod.iter()
            .filter(|h| !h.is_reference && h.bases.as_slice() == BASES_A)
            .count(),
        1
    );
}

/// Hamming-1 SNPs are distinct sequences and must both survive.
#[test]
fn forensic_6r81_distinct_sequences_are_not_collapsed() {
    let loc = GenomeLoc::new(1, BASES_A.len() as u64);
    assert_eq!(BASES_A.len(), BASES_B.len());
    let diffs: Vec<usize> = BASES_A
        .iter()
        .zip(BASES_B.iter())
        .enumerate()
        .filter(|(_, (a, b))| a != b)
        .map(|(i, _)| i)
        .collect();
    assert_eq!(diffs, vec![19]);
    let haps = vec![
        hap_match(BASES_A, true, loc),
        hap_match(BASES_B, false, loc),
    ];
    let prod = production_trim_to(haps.clone(), BASES_A, 1, BASES_A.len() as u64);
    let java = java_trim_down_haplotypes(&haps, &loc);
    assert_eq!(java.len(), 2);
    assert_eq!(prod.len(), 2);
    assert_eq!(unique_seqs(&prod).len(), 2);
    assert!(prod
        .iter()
        .any(|h| h.is_reference && h.bases.as_slice() == BASES_A));
    assert!(prod
        .iter()
        .any(|h| !h.is_reference && h.bases.as_slice() == BASES_B));
}

/// Live kernel: unique-sequence SET already matched Java; extra column was the
/// ref+alt duplicate. After 6R.81, object count must match unique count.
#[test]
fn live_pairhmm_haplotype_population_matches_java_object_count() {
    if std::env::var("HOLDOUT_6R81").ok().as_deref() != Some("1") {
        eprintln!("skip: set HOLDOUT_6R81=1");
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

    let rust_dump_path = std::env::temp_dir().join("6r81_rust_pairhmm_inputs.txt");
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
    let legacy = legacy_trim_keep_ref_state(&untrimmed.haplotypes, &span);

    let mut legacy_dupes: Vec<(String, usize, usize)> = Vec::new();
    {
        let mut by_bases: HashMap<Vec<u8>, Vec<bool>> = HashMap::new();
        for h in &legacy {
            by_bases
                .entry(h.bases.clone())
                .or_default()
                .push(h.is_reference);
        }
        for (bases, flags) in by_bases {
            if flags.len() > 1 {
                legacy_dupes.push((
                    String::from_utf8_lossy(&bases).into_owned(),
                    bases.len(),
                    flags.len(),
                ));
            }
        }
    }
    legacy_dupes.sort_by(|a, b| a.0.cmp(&b.0));

    eprintln!(
        "6R.81 untrimmed n={} unique={}",
        untrimmed.haplotypes.len(),
        unique_seqs(&untrimmed.haplotypes).len()
    );
    eprintln!(
        "6R.81 trim objects: legacy={} java_replay={} production={} unique_prod={}",
        legacy.len(),
        java_replay.len(),
        rust_trim.haplotypes.len(),
        unique_seqs(&rust_trim.haplotypes).len()
    );
    eprintln!(
        "6R.81 legacy same-bases object groups={}",
        legacy_dupes.len()
    );
    for (seq, len, n) in &legacy_dupes {
        eprintln!(
            "  dupe n={} len={} fnv={:016x} seq={}",
            n,
            len,
            fnv1a64(seq.as_bytes()),
            seq
        );
    }

    assert_eq!(
        unique_seqs(&legacy).len(),
        unique_seqs(&java_replay).len(),
        "unique-sequence SET already matched; extra is object duplication"
    );
    assert_eq!(
        legacy.len(),
        java_replay.len() + legacy_dupes.iter().map(|d| d.2 - 1).sum::<usize>(),
        "legacy extra objects are exactly the same-bases duplicates"
    );
    assert_eq!(
        rust_trim.haplotypes.len(),
        java_replay.len(),
        "production trim_to object count must match Java trimDown"
    );
    assert_eq!(
        unique_seqs(&rust_trim.haplotypes),
        unique_seqs(&java_replay)
    );

    if rust_dump_path.is_file() && java_dump_path.is_file() {
        let rust_dump = std::fs::read_to_string(&rust_dump_path).expect("rust dump");
        let java_dump = std::fs::read_to_string(&java_dump_path).expect("java dump");
        let rust_haps = snp_motif_haps(&rust_dump, MOTIF);
        let java_haps = snp_motif_haps(&java_dump, MOTIF);
        let rs: HashSet<&str> = rust_haps.iter().map(|s| s.as_str()).collect();
        let js: HashSet<&str> = java_haps.iter().map(|s| s.as_str()).collect();
        let only_j = js.difference(&rs).count();
        let only_r = rs.difference(&js).count();
        eprintln!(
            "6R.81 PairHMM columns java={} rust={} unique_java={} unique_rust={} JAVA_ONLY={} RUST_ONLY={}",
            java_haps.len(),
            rust_haps.len(),
            js.len(),
            rs.len(),
            only_j,
            only_r
        );
        assert_eq!(only_j, 0, "JAVA_ONLY sequences");
        assert_eq!(only_r, 0, "RUST_ONLY sequences");
        assert_eq!(java_haps.len(), rust_haps.len(), "PairHMM column count");
        assert_eq!(
            java_haps.len(),
            js.len(),
            "Java columns are unique sequences"
        );
        assert_eq!(
            rust_haps.len(),
            rs.len(),
            "Rust columns are unique sequences"
        );
    }
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
