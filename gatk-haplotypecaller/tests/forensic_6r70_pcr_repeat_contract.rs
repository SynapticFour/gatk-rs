//! 6R.70 coordinate-free: Java PCR cache vs `findTandemRepeatUnits`.
//!
//! Java 4.4.0.0 (`2dbc025821bc5f686c423ff332a41e6cef892a77`)
//! `ReadLikelihoodCalculationEngine.findTandemRepeatUnits` +
//! `PairHMMLikelihoodCalculationEngine.applyPCRErrorModel`.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r70_pcr_repeat_contract
//! ```

use gatk_haplotypecaller::pcr_error_model::{error_model_adjusted_qual, PcrErrorModel};

const JAVA_MAX_STR_UNIT: usize = 8;
const JAVA_MAX_REPEAT: usize = 20;

/// GATK 4.4 `MathUtils.fastRound`: `(int)(d + 0.5)` for positive `d`.
fn java_fast_round(d: f64) -> i32 {
    if d > 0.0 {
        (d + 0.5) as i32
    } else {
        (d - 0.5) as i32
    }
}

/// GATK 4.4 CONSERVATIVE `getErrorModelAdjustedQual`.
fn java_44_cache(r: usize) -> u8 {
    let q = 40.0 - (r as f64 / (3.0 * std::f64::consts::PI)).exp() + 1.0;
    java_fast_round(q).max(10) as u8
}

fn rust_conservative_cache(r: usize) -> u8 {
    error_model_adjusted_qual(r, PcrErrorModel::Conservative.rate_factor().unwrap())
}

fn equal_range(a: &[u8], a_off: usize, b: &[u8], b_off: usize, len: usize) -> bool {
    a.get(a_off..a_off + len) == b.get(b_off..b_off + len)
}

/// GATK 4.4 `GATKVariantContextUtils.findNumberOfRepetitions` (subarray form).
fn find_number_of_repetitions(
    unit: &[u8],
    unit_off: usize,
    unit_len: usize,
    test: &[u8],
    test_off: usize,
    test_len: usize,
    leading: bool,
) -> usize {
    if unit_len == 0 || test_len == 0 {
        return 0;
    }
    let length_difference = test_len as isize - unit_len as isize;
    if leading {
        let mut n = 0usize;
        let mut start = 0isize;
        while start <= length_difference {
            if equal_range(test, (start as usize) + test_off, unit, unit_off, unit_len) {
                n += 1;
                start += unit_len as isize;
            } else {
                return n;
            }
        }
        n
    } else {
        let mut n = 0usize;
        let mut start = length_difference;
        while start >= 0 {
            if equal_range(test, (start as usize) + test_off, unit, unit_off, unit_len) {
                n += 1;
                start -= unit_len as isize;
            } else {
                return n;
            }
        }
        n
    }
}

/// GATK 4.4 `ReadLikelihoodCalculationEngine.findTandemRepeatUnits`.
/// Returns (repeat unit, repeat length used as PCR cache index).
fn java_find_tandem_repeat_units(read: &[u8], offset: usize) -> (Vec<u8>, usize) {
    assert!(offset < read.len());
    let mut max_bw = 0usize;
    let mut best_bw: Vec<u8> = vec![read[offset]];
    for str_len in 1..=JAVA_MAX_STR_UNIT {
        if (offset + 1).checked_sub(str_len).is_none() {
            break;
        }
        max_bw = find_number_of_repetitions(
            read,
            offset + 1 - str_len,
            str_len,
            read,
            0,
            offset + 1,
            false,
        );
        if max_bw > 1 {
            best_bw = read[offset + 1 - str_len..=offset].to_vec();
            break;
        }
    }
    let mut best_unit = best_bw.clone();
    let mut max_rl = max_bw;

    if offset < read.len() - 1 {
        let mut best_fw: Vec<u8> = vec![read[offset + 1]];
        let mut max_fw = 0usize;
        for str_len in 1..=JAVA_MAX_STR_UNIT {
            if offset + str_len + 1 > read.len() {
                break;
            }
            max_fw = find_number_of_repetitions(
                read,
                offset + 1,
                str_len,
                read,
                offset + 1,
                read.len() - offset - 1,
                true,
            );
            if max_fw > 1 {
                best_fw = read[offset + 1..offset + 1 + str_len].to_vec();
                break;
            }
        }
        if best_fw == best_bw {
            max_rl = max_bw + max_fw;
            best_unit = best_fw;
        } else {
            let test = &read[..=offset];
            max_bw =
                find_number_of_repetitions(&best_fw, 0, best_fw.len(), test, 0, test.len(), false);
            max_rl = max_fw + max_bw;
            best_unit = best_fw;
        }
    }
    if max_rl > JAVA_MAX_REPEAT {
        max_rl = JAVA_MAX_REPEAT;
    }
    (best_unit, max_rl)
}

/// 6R.70 historical Rust: homopolymer-only walk. Production is the Java STR
/// finder as of 6R.71; this helper keeps the 6R.70 divergence table intact.
fn homopolymer_repeat_units(bases: &[u8], i: usize) -> usize {
    if bases.is_empty() || i >= bases.len() {
        return 0;
    }
    let base = bases[i];
    let mut len = 1usize;
    let mut j = i;
    while j > 0 && bases[j - 1] == base {
        len += 1;
        j -= 1;
    }
    j = i;
    while j + 1 < bases.len() && bases[j + 1] == base {
        len += 1;
        j += 1;
    }
    len
}

fn rust_unit_at(read: &[u8], offset: usize) -> (Vec<u8>, usize) {
    if read.is_empty() || offset >= read.len() {
        return (Vec::new(), 0);
    }
    (vec![read[offset]], homopolymer_repeat_units(read, offset))
}

/// Pinned Java CONSERVATIVE cache for r=0..=20. Index 20 is the finder clip;
/// `applyPCRErrorModel` then uses `cache[repeatLength]` with no extra min.
const JAVA_CACHE_0_TO_20: [u8; 21] = [
    40, 40, 40, 40, 39, 39, 39, 39, 39, 38, 38, 38, 37, 37, 37, 36, 36, 35, 34, 33, 33,
];

#[test]
fn java_conservative_cache_matches_rust_at_every_index_0_to_20() {
    for r in 0..=20 {
        let j = java_44_cache(r);
        let rs = rust_conservative_cache(r);
        assert_eq!(j, JAVA_CACHE_0_TO_20[r], "Java cache[{r}]");
        assert_eq!(rs, j, "r={r} Java {j} Rust {rs}");
        eprintln!("cache r={r} Java={j} Rust={rs} equal=yes");
    }
    assert_eq!(java_44_cache(0), 40);
    assert_eq!(java_44_cache(1), 40);
    assert_eq!(rust_conservative_cache(0), 40);
    assert_eq!(rust_conservative_cache(1), 40);
    // Java min-10 is first hit at r=33, past the finder clip of 20.
    assert_eq!(java_44_cache(32), 11);
    assert_eq!(java_44_cache(33), 10);
    for r in 33..=40 {
        assert_eq!(java_44_cache(r), 10);
        assert_eq!(rust_conservative_cache(r), 10);
    }
    assert_eq!(rust_conservative_cache(20), 33);
}

#[test]
fn find_number_of_repetitions_matches_java_unit_examples() {
    // GATKVariantContextUtilsUnitTest vectors at the 4.4 pin.
    assert_eq!(
        find_number_of_repetitions(b"AT", 0, 2, b"GATAT", 0, 5, false),
        2
    );
    assert_eq!(
        find_number_of_repetitions(b"AT", 0, 2, b"GATAT", 0, 5, true),
        0
    );
    assert_eq!(
        find_number_of_repetitions(b"AT", 0, 2, b"ATATG", 0, 5, true),
        2
    );
    assert_eq!(find_number_of_repetitions(b"T", 0, 1, b"T", 0, 1, true), 1);
}

#[test]
fn rust_round_vs_java_fast_round_do_not_explain_the_cache_gap() {
    // Same formula body, only INITIAL/rate/min differ. Rounding is not the fork.
    for r in 0..=20 {
        let q40 = 40.0 - (r as f64 / (3.0 * std::f64::consts::PI)).exp() + 1.0;
        let fast = java_fast_round(q40);
        let rnd = q40.round() as i32;
        assert_eq!(fast, rnd, "r={r} fastRound vs f64::round on Java formula");
    }
}

#[test]
fn homopolymer_repeat_lengths_agree() {
    for seq in [b"AAAAAA".as_slice(), b"TTTTTT".as_slice()] {
        for pos in 0..seq.len() {
            let (ju, jl) = java_find_tandem_repeat_units(seq, pos);
            let (ru, rl) = rust_unit_at(seq, pos);
            assert_eq!(
                jl,
                rl,
                "{} pos {pos}: Java {jl} Rust {rl}",
                std::str::from_utf8(seq).unwrap()
            );
            assert_eq!(ju, ru);
        }
    }
}

#[test]
fn dinucleotide_and_trinucleotide_java_is_broader_than_homopolymer() {
    let at = b"ATATAT";
    let (ju0, jl0) = java_find_tandem_repeat_units(at, 0);
    let (_, rl0) = rust_unit_at(at, 0);
    assert_eq!(ju0, b"TA");
    assert_eq!(jl0, 2);
    assert_eq!(rl0, 1, "Rust homopolymer at ATATAT[0] is 1");

    let ac = b"ACACACAC";
    let (ju0, jl0) = java_find_tandem_repeat_units(ac, 0);
    let (_, rl0) = rust_unit_at(ac, 0);
    assert!(
        jl0 >= 2,
        "Java ACAC… at 0 should be multi-copy AC/CA, got {jl0}"
    );
    assert_eq!(ju0.len(), 2);
    assert_eq!(rl0, 1);

    let acg = b"ACGACGACG";
    let (ju0, jl0) = java_find_tandem_repeat_units(acg, 0);
    let (_, rl0) = rust_unit_at(acg, 0);
    assert_eq!(ju0, b"CGA");
    assert_eq!(jl0, 2);
    assert_eq!(rl0, 1);
}

#[test]
fn non_repeat_and_boundaries() {
    let seq = b"ACGT";
    for pos in 0..seq.len() {
        let (_, jl) = java_find_tandem_repeat_units(seq, pos);
        let (_, rl) = rust_unit_at(seq, pos);
        assert_eq!(rl, 1);
        // Java leftover maxFW from the last STR width can be 0; GOP cache[0] and cache[1] are both Q40.
        assert!(
            jl <= 1,
            "ACGT pos {pos}: Java len {jl} should not invent a tandem"
        );
    }

    let one = b"T";
    let (ju, jl) = java_find_tandem_repeat_units(one, 0);
    let (ru, rl) = rust_unit_at(one, 0);
    assert_eq!(ju, b"T");
    assert_eq!(ru, b"T");
    assert_eq!(jl, 1);
    assert_eq!(rl, 1);

    let interrupted = b"AAAATAAAA";
    // Last A of the left run: Rust counts the homopolymer (4); Java switches to the
    // forward unit T and reports length 1.
    let (ju3, j3) = java_find_tandem_repeat_units(interrupted, 3);
    let (_, r3) = rust_unit_at(interrupted, 3);
    assert_eq!(ju3, b"T");
    assert_eq!(j3, 1);
    assert_eq!(r3, 4);
    // Interrupting T: Rust homopolymer 1; Java takes the forward A-run of 4.
    let (ju4, j4) = java_find_tandem_repeat_units(interrupted, 4);
    let (_, r4) = rust_unit_at(interrupted, 4);
    assert_eq!(ju4, b"A");
    assert_eq!(j4, 4);
    assert_eq!(r4, 1, "T interrupting A-runs is a homopolymer of 1");
}

fn row(seq: &[u8], pos: usize) -> String {
    let (ju, jl) = java_find_tandem_repeat_units(seq, pos);
    let (ru, rl) = rust_unit_at(seq, pos);
    format!(
        "{} | {pos} | {} | {jl} | {} | {rl}",
        std::str::from_utf8(seq).unwrap(),
        std::str::from_utf8(&ju).unwrap(),
        std::str::from_utf8(&ru).unwrap()
    )
}

#[test]
fn repeat_fixture_table_java_vs_rust() {
    eprintln!("sequence | pos | Java unit | Java len | Rust unit | Rust len");
    let fixtures: &[&[u8]] = &[
        b"AAAAAA",
        b"TTTTTT",
        b"ATATAT",
        b"ACACACAC",
        b"ACGACGACG",
        b"ACGT",
        b"T",
        b"AAAATAAAA",
        b"TGCA",
        b"TAAGAAAA", // canonical-class prefix (read 0 / base 0)
    ];
    for seq in fixtures {
        for pos in 0..seq.len() {
            eprintln!("{}", row(seq, pos));
            let (_, jl) = java_find_tandem_repeat_units(seq, pos);
            assert!(
                jl >= 1,
                "cache[0] not selected by the finder on nonempty reads"
            );
            assert!(jl <= JAVA_MAX_REPEAT);
        }
    }

    // Homopolymers: lengths agree (Case A on length). Cache now matches (6R.72).
    for pos in 0..6 {
        assert_eq!(java_find_tandem_repeat_units(b"AAAAAA", pos).1, 6);
        assert_eq!(homopolymer_repeat_units(b"AAAAAA", pos), 6);
    }

    // Multi-base tandem: Java is broader (Case B).
    assert_eq!(
        java_find_tandem_repeat_units(b"ATATAT", 0),
        (b"TA".to_vec(), 2)
    );
    assert_eq!(homopolymer_repeat_units(b"ATATAT", 0), 1);
    assert_eq!(
        java_find_tandem_repeat_units(b"ATATAT", 1),
        (b"AT".to_vec(), 3)
    );
    assert_eq!(homopolymer_repeat_units(b"ATATAT", 1), 1);

    // Single-base and start-of-read T: Java index is 1, not 0.
    assert_eq!(java_find_tandem_repeat_units(b"T", 0), (b"T".to_vec(), 1));
    assert_eq!(
        java_find_tandem_repeat_units(b"TGCA", 0),
        (b"G".to_vec(), 1)
    );
    assert_eq!(homopolymer_repeat_units(b"TGCA", 0), 1);

    // Canonical-class prefix: FW A-run of 2, not homopolymer T.
    let canon = b"TAAGAAAA";
    let (ju, jl) = java_find_tandem_repeat_units(canon, 0);
    assert_eq!(ju, b"A");
    assert_eq!(jl, 2);
    assert_eq!(homopolymer_repeat_units(canon, 0), 1);
    assert_eq!(java_44_cache(1), 40);
    assert_eq!(java_44_cache(2), 40);
    assert_eq!(44u8.min(java_44_cache(2)), 40);
    assert_eq!(44u8.min(rust_conservative_cache(1)), 40);
    // Length 1 vs 2 stays GOP-silent: cache[1]==cache[2]==40. Finder proof is unit/length.

    // Last base is findable but never written by applyPCRErrorModel.
    let last = b"ACGT".len() - 1;
    assert_eq!(java_find_tandem_repeat_units(b"ACGT", last).1, 1);
}

#[test]
fn last_base_is_never_pcr_indexed_by_the_apply_loop() {
    // Java `for (i = 1; i < len; i++)` uses offset i-1. Last index is not a cache lookup.
    let seq = b"ACGT";
    let last = seq.len() - 1;
    assert_eq!(last, 3);
    let (_, jl) = java_find_tandem_repeat_units(seq, last);
    let (_, rl) = rust_unit_at(seq, last);
    assert_eq!(rl, 1);
    let _ = jl;
}

#[test]
fn cache_index_1_and_0_are_both_java_q40() {
    assert_eq!(java_44_cache(0), 40);
    assert_eq!(java_44_cache(1), 40);
    assert_eq!(java_44_cache(2), 40);
    assert_eq!(java_44_cache(3), 40);
}

#[test]
fn after_matching_repeat_length_one_gop_matches_via_cache() {
    // Isolated cache proof when lengths *do* agree (homopolymer / unique).
    assert_eq!(44u8.min(java_44_cache(1)), 40);
    assert_eq!(44u8.min(rust_conservative_cache(1)), 40);
    assert_eq!(44u8.min(java_44_cache(2)), 40);
    assert_eq!(44u8.min(rust_conservative_cache(2)), 40);
}

#[test]
fn pcr_loop_uses_offset_i_minus_one_before_bq_cap() {
    // Documented Java order: PCR min-cap, then BQ/MAPQ cap, then IQ/DQ floor, GCP later.
    let seq = b"ACGTACGT";
    let mut java_gop = vec![44u8; seq.len()];
    for i in 1..seq.len() {
        let (_, rl) = java_find_tandem_repeat_units(seq, i - 1);
        java_gop[i - 1] = java_gop[i - 1].min(java_44_cache(rl.min(JAVA_MAX_REPEAT)));
    }
    assert_eq!(java_gop[seq.len() - 1], 44, "last base not PCR-written");
    assert_eq!(java_gop[0], 40);
}
