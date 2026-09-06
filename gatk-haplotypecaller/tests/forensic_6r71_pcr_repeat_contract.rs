//! 6R.71 coordinate-free: Java 4.4 `findTandemRepeatUnits` vs production Rust.
//!
//! SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`.
//! Repeat *length* is the unit *count* (not span in bases), clipped to 20 after
//! combining backward+forward. Candidate widths are tried **1..=8, first count>1**.
//!
//! Finder is frozen. 6R.72 PCR cache matches Java; canonical GOP is 40 after
//! BI-sub. Production fill Q45 is still a separate 6R.68 fork.
//!
//! ```text
//! cargo test -p gatk-haplotypecaller --test forensic_6r71_pcr_repeat_contract
//! ```

use gatk_haplotypecaller::pcr_error_model::{
    error_model_adjusted_qual, find_tandem_repeat_units, tandem_repeat_units, PcrErrorModel,
    MAX_REPEAT_LENGTH, MAX_STR_UNIT_LENGTH,
};

const JAVA_MAX_STR_UNIT: usize = 8;
const JAVA_MAX_REPEAT: usize = 20;

fn java_fast_round(d: f64) -> i32 {
    if d > 0.0 {
        (d + 0.5) as i32
    } else {
        (d - 0.5) as i32
    }
}

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
            if equal_range(test, start as usize + test_off, unit, unit_off, unit_len) {
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
            if equal_range(test, start as usize + test_off, unit, unit_off, unit_len) {
                n += 1;
                start -= unit_len as isize;
            } else {
                return n;
            }
        }
        n
    }
}

/// Pinned Java `ReadLikelihoodCalculationEngine.findTandemRepeatUnits`.
fn java_find_tandem_repeat_units(read: &[u8], offset: usize) -> (Vec<u8>, usize) {
    assert!(offset < read.len());
    let mut max_bw = 0usize;
    let mut best_bw: Vec<u8> = vec![read[offset]];
    for str_len in 1..=JAVA_MAX_STR_UNIT {
        if offset + 1 < str_len {
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

fn assert_matches_java(seq: &[u8], pos: usize) {
    let want = java_find_tandem_repeat_units(seq, pos);
    let got = find_tandem_repeat_units(seq, pos);
    assert_eq!(
        got,
        want,
        "{} pos {pos}: production {:?} Java {:?}",
        std::str::from_utf8(seq).unwrap_or("?"),
        (std::str::from_utf8(&got.0).unwrap_or("?"), got.1),
        (std::str::from_utf8(&want.0).unwrap_or("?"), want.1)
    );
    assert_eq!(tandem_repeat_units(seq, pos), want.1);
}

#[test]
fn constants_match_pinned_java() {
    assert_eq!(MAX_STR_UNIT_LENGTH, 8);
    assert_eq!(MAX_REPEAT_LENGTH, 20);
    assert_eq!(JAVA_MAX_STR_UNIT, 8);
    assert_eq!(JAVA_MAX_REPEAT, 20);
}

#[test]
fn width_1_homopolymers() {
    for seq in [b"AAAAAA".as_slice(), b"TTTTTT".as_slice()] {
        for pos in 0..seq.len() {
            assert_matches_java(seq, pos);
            let (u, n) = find_tandem_repeat_units(seq, pos);
            assert_eq!(u, vec![seq[pos]]);
            assert_eq!(n, 6);
        }
    }
}

#[test]
fn width_gt1_observed_java_units_not_generic_str() {
    // Do not assume ATATAT → AT/3. Java at [0] is TA/2 (FW unit, leftover T dropped).
    assert_eq!(
        java_find_tandem_repeat_units(b"ATATAT", 0),
        (b"TA".to_vec(), 2)
    );
    assert_eq!(find_tandem_repeat_units(b"ATATAT", 0), (b"TA".to_vec(), 2));
    assert_eq!(
        java_find_tandem_repeat_units(b"ATATAT", 1),
        (b"AT".to_vec(), 3)
    );
    assert_eq!(find_tandem_repeat_units(b"ATATAT", 1), (b"AT".to_vec(), 3));
    for pos in 0..6 {
        assert_matches_java(b"ATATAT", pos);
    }

    assert_eq!(
        java_find_tandem_repeat_units(b"ACACACAC", 0),
        (b"CA".to_vec(), 3)
    );
    for pos in 0..8 {
        assert_matches_java(b"ACACACAC", pos);
    }

    assert_eq!(
        java_find_tandem_repeat_units(b"ACGACGACG", 0),
        (b"CGA".to_vec(), 2)
    );
    for pos in 0..9 {
        assert_matches_java(b"ACGACGACG", pos);
    }
}

#[test]
fn aaaataaaa_java_is_not_longest_homopolymer() {
    let seq = b"AAAATAAAA";
    // [3] last A of the left run: FW unit T does not match BW A, so FW wins → T/1.
    assert_eq!(java_find_tandem_repeat_units(seq, 3), (b"T".to_vec(), 1));
    assert_eq!(find_tandem_repeat_units(seq, 3), (b"T".to_vec(), 1));
    // [4] interrupting T: FW A-run of 4.
    assert_eq!(java_find_tandem_repeat_units(seq, 4), (b"A".to_vec(), 4));
    assert_eq!(find_tandem_repeat_units(seq, 4), (b"A".to_vec(), 4));
    for pos in 0..seq.len() {
        assert_matches_java(seq, pos);
    }
}

#[test]
fn canonical_prefix_anti_masking_unit_and_length() {
    let seq = b"TAAGAAAA";
    let (ju, jl) = java_find_tandem_repeat_units(seq, 0);
    let (ru, rl) = find_tandem_repeat_units(seq, 0);
    assert_eq!(ju, b"A");
    assert_eq!(jl, 2);
    assert_eq!(ru, b"A");
    assert_eq!(rl, 2);
    // GOP cannot prove the finder: Java cache[1]==cache[2]==40 (and Rust now too).
    assert_eq!(java_44_cache(1), 40);
    assert_eq!(java_44_cache(2), 40);
    assert_eq!(44u8.min(java_44_cache(2)), 40);
    assert_eq!(44u8.min(rust_conservative_cache(2)), 40);
    assert_eq!(rust_conservative_cache(2), java_44_cache(2));
}

#[test]
fn first_width_with_count_gt_one_wins() {
    // Homopolymer: width 1 matches first even though AA also tiles.
    assert_eq!(find_tandem_repeat_units(b"AAAAAA", 2), (b"A".to_vec(), 6));

    // Width 3 matches before width 6 on a 3-mer tandem.
    let seq = b"ACGACGACGACG";
    let (u, n) = find_tandem_repeat_units(seq, 0);
    assert_eq!(u, b"CGA");
    assert_eq!(n, 3);
    assert_ne!(u.as_slice(), b"CGACGA");
    assert_matches_java(seq, 0);

    // ATATAT[0]: width 1 FW count is 1, so width 2 (TA) is the first >1.
    assert_eq!(find_tandem_repeat_units(b"ATATAT", 0).0, b"TA");
}

#[test]
fn clip_is_total_unit_count_not_a_20_base_search_window() {
    // 25-base homopolymer: finder walks the full sequence, then clips the *count* to 20.
    let long_a = vec![b'A'; 25];
    let (u0, n0) = find_tandem_repeat_units(&long_a, 0);
    assert_eq!(u0, b"A");
    assert_eq!(n0, 20, "pos 0: 1+24 units clipped to 20");
    assert_eq!(java_find_tandem_repeat_units(&long_a, 0).1, 20);

    let (um, nm) = find_tandem_repeat_units(&long_a, 12);
    assert_eq!(um, b"A");
    assert_eq!(nm, 20, "middle: 13+12=25 clipped to 20");

    let last = long_a.len() - 1;
    let (ul, nl) = find_tandem_repeat_units(&long_a, last);
    assert_eq!(ul, b"A");
    assert_eq!(
        nl, 20,
        "last base: BW 25 clipped to 20 (finder, not PCR loop)"
    );

    // 21 copies of AT (42 bases): unit count 21 → clip 20, which is 40 bases of span.
    let mut at = Vec::new();
    for _ in 0..21 {
        at.extend_from_slice(b"AT");
    }
    let (u, n) = find_tandem_repeat_units(&at, 1);
    assert_eq!(n, 20);
    assert_eq!(u.len(), 2);
    assert_eq!(java_find_tandem_repeat_units(&at, 1).1, 20);
}

#[test]
fn boundaries_and_no_repeat() {
    let seq = b"ACGT";
    for pos in 0..seq.len() {
        assert_matches_java(seq, pos);
        assert_eq!(find_tandem_repeat_units(seq, pos).1, 1);
    }
    assert_eq!(find_tandem_repeat_units(b"T", 0), (b"T".to_vec(), 1));

    let hom = b"AAAA";
    assert_matches_java(hom, 0);
    assert_matches_java(hom, 1);
    assert_matches_java(hom, hom.len() - 2);
    assert_matches_java(hom, hom.len() - 1);
}

#[test]
fn last_base_rule_stays_in_the_pcr_loop_not_the_finder() {
    let seq = b"AAAAAA";
    let last = seq.len() - 1;
    assert_eq!(find_tandem_repeat_units(seq, last).1, 6);
    // applyPCRErrorModel does not write ins[last]; that is the caller, not this function.
}

#[test]
fn different_java_cache_slots_make_wrong_length_observable() {
    // Java cache[1]=40, cache[4]=39. AAAATAAAA[4] is length 4.
    assert_eq!(java_44_cache(1), 40);
    assert_eq!(java_44_cache(4), 39);
    let (_, n) = find_tandem_repeat_units(b"AAAATAAAA", 4);
    assert_eq!(n, 4);
    assert_eq!(44u8.min(java_44_cache(n)), 39);
    assert_eq!(44u8.min(java_44_cache(1)), 40);
    assert_eq!(rust_conservative_cache(1), 40);
    assert_eq!(rust_conservative_cache(4), 39);
}

#[test]
fn fixture_table_production_equals_java() {
    let fixtures: &[&[u8]] = &[
        b"AAAAAA",
        b"TTTTTT",
        b"ATATAT",
        b"ACACACAC",
        b"ACGACGACG",
        b"AAAATAAAA",
        b"TAAGAAAA",
        b"ACGT",
        b"T",
        b"TGCA",
    ];
    eprintln!("sequence | pos | Java unit | Java len | Rust unit | Rust len");
    for seq in fixtures {
        for pos in 0..seq.len() {
            let (ju, jl) = java_find_tandem_repeat_units(seq, pos);
            let (ru, rl) = find_tandem_repeat_units(seq, pos);
            eprintln!(
                "{} | {pos} | {} | {jl} | {} | {rl}",
                std::str::from_utf8(seq).unwrap(),
                std::str::from_utf8(&ju).unwrap(),
                std::str::from_utf8(&ru).unwrap()
            );
            assert_eq!((ru, rl), (ju, jl));
        }
    }
}
