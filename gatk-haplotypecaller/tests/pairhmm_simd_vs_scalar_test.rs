//! SIMD / packed Logless PairHMM vs scalar Logless (and Log10 within policy).

use gatk_haplotypecaller::{
    best_simd_backend, log10_pairhmm_likelihood_exact, logless_pairhmm_likelihood, pairhmm_fp_eq,
    resolve_pair_hmm_impl, score_read_haps_logless, PairHmmFpPolicy, PairHmmImpl,
    GATK_PARITY_DEFAULT_DEL_QUAL, GATK_PARITY_DEFAULT_GCP, GATK_PARITY_DEFAULT_INS_QUAL,
};

fn parity_quals(n: usize) -> (Vec<u8>, Vec<u8>, Vec<u8>, Vec<u8>) {
    (
        vec![30u8; n],
        vec![GATK_PARITY_DEFAULT_INS_QUAL; n],
        vec![GATK_PARITY_DEFAULT_DEL_QUAL; n],
        vec![GATK_PARITY_DEFAULT_GCP; n],
    )
}

fn bases_pattern(len: usize, seed: u64) -> Vec<u8> {
    const ALPHA: &[u8] = b"ACGT";
    (0..len)
        .map(|i| ALPHA[((seed.wrapping_mul(1103515245).wrapping_add(i as u64)) % 4) as usize])
        .collect()
}

#[test]
fn simd_matches_scalar_logless_on_synthetic_suite() {
    let backend = resolve_pair_hmm_impl(PairHmmImpl::Simd);
    let policy = PairHmmFpPolicy {
        abs_epsilon: 1e-9,
        rel_epsilon: 1e-8,
    };
    let lengths = [1usize, 2, 3, 7, 8, 15, 16, 31, 32, 64, 100, 151, 300];
    let mut cases = 0usize;
    for &rn in &lengths {
        for &hn in &lengths {
            for seed in 0..3u64 {
                let read = bases_pattern(rn, seed);
                let (quals, ins, del, gcp) = parity_quals(rn);
                let mut haps: Vec<Vec<u8>> = Vec::new();
                // Mix equal-length and uneven packs (SIMD width edges).
                for k in 0..5u64 {
                    let mut h = bases_pattern(hn, seed.wrapping_add(k * 17));
                    if k % 2 == 1 && hn > 1 {
                        h.truncate(hn - 1);
                    }
                    if k % 3 == 0 && !h.is_empty() {
                        h[0] = b'T';
                    }
                    haps.push(h);
                }
                let hap_refs: Vec<&[u8]> = haps.iter().map(|h| h.as_slice()).collect();
                let simd =
                    score_read_haps_logless(backend, &read, &quals, &hap_refs, &ins, &del, &gcp)
                        .expect("simd");
                for (i, hap) in hap_refs.iter().enumerate() {
                    let scalar = logless_pairhmm_likelihood(&read, &quals, hap, &ins, &del, &gcp)
                        .expect("scalar");
                    assert!(
                        pairhmm_fp_eq(simd[i], scalar, policy)
                            || (simd[i].is_infinite() && scalar.is_infinite()),
                        "mismatch rn={rn} hn={} seed={seed} i={i}: simd={} scalar={} backend={backend:?}",
                        hap.len(),
                        simd[i],
                        scalar
                    );
                    cases += 1;
                }
            }
        }
    }
    assert!(cases >= 200, "expected hundreds of cases, got {cases}");
    let _ = best_simd_backend();
}

/// Phenotype-class differential: read lengths 100–300 and hap packs 8–128
/// (showcase Phase B; not locus pins).
#[test]
fn simd_matches_scalar_on_read_len_hap_count_phenotypes() {
    let backend = resolve_pair_hmm_impl(PairHmmImpl::Simd);
    let policy = PairHmmFpPolicy {
        abs_epsilon: 1e-9,
        rel_epsilon: 1e-8,
    };
    let mut cases = 0usize;
    for &read_len in &[100usize, 200, 300] {
        for &hap_count in &[8usize, 32, 64, 128] {
            let read = bases_pattern(read_len, 7);
            let (quals, ins, del, gcp) = parity_quals(read_len);
            let haps: Vec<Vec<u8>> = (0..hap_count)
                .map(|k| {
                    let mut h = bases_pattern(read_len + (k % 5), 11 + k as u64);
                    if k % 4 == 0 && !h.is_empty() {
                        h[0] = b'G';
                    }
                    h
                })
                .collect();
            let hap_refs: Vec<&[u8]> = haps.iter().map(|h| h.as_slice()).collect();
            let simd = score_read_haps_logless(backend, &read, &quals, &hap_refs, &ins, &del, &gcp)
                .expect("simd phenotype");
            assert_eq!(simd.len(), hap_count);
            for (i, hap) in hap_refs.iter().enumerate() {
                let scalar = logless_pairhmm_likelihood(&read, &quals, hap, &ins, &del, &gcp)
                    .expect("scalar phenotype");
                assert!(
                    pairhmm_fp_eq(simd[i], scalar, policy)
                        || (simd[i].is_infinite() && scalar.is_infinite()),
                    "phenotype rn={read_len} haps={hap_count} i={i}: simd={} scalar={}",
                    simd[i],
                    scalar
                );
                cases += 1;
            }
        }
    }
    assert!(cases >= 8 + 32 + 64 + 128, "got {cases}");
}

#[test]
fn logless_and_log10_both_finite_same_rank_order() {
    // Logless vs Log10 are different engines (Java likewise). Absolute Δ can be
    // O(0.1–1) log10 on short synthetic pairs; we only require finiteness and that
    // a perfect match outscores a mismatch for both engines.
    let read = b"ACGTACGTAC";
    let (quals, ins, del, gcp) = parity_quals(read.len());
    let hap_ok = b"ACGTACGTAC";
    let hap_bad = b"ACGTACGTAT";
    let ll_ok = logless_pairhmm_likelihood(read, &quals, hap_ok, &ins, &del, &gcp).unwrap();
    let ll_bad = logless_pairhmm_likelihood(read, &quals, hap_bad, &ins, &del, &gcp).unwrap();
    let l10_ok = log10_pairhmm_likelihood_exact(read, &quals, hap_ok, &ins, &del, &gcp).unwrap();
    let l10_bad = log10_pairhmm_likelihood_exact(read, &quals, hap_bad, &ins, &del, &gcp).unwrap();
    assert!(ll_ok.is_finite() && ll_bad.is_finite());
    assert!(l10_ok.is_finite() && l10_bad.is_finite());
    assert!(ll_ok > ll_bad, "logless rank: {ll_ok} vs {ll_bad}");
    assert!(l10_ok > l10_bad, "log10 rank: {l10_ok} vs {l10_bad}");
}

#[test]
fn f32_retry_path_matches_scalar_within_looser_policy() {
    let backend = resolve_pair_hmm_impl(PairHmmImpl::SimdF32);
    let policy = PairHmmFpPolicy {
        abs_epsilon: 1e-3,
        rel_epsilon: 1e-2,
    };
    let read = bases_pattern(80, 42);
    let (quals, ins, del, gcp) = parity_quals(read.len());
    let haps: Vec<Vec<u8>> = (0..8).map(|k| bases_pattern(90, 100 + k)).collect();
    let hap_refs: Vec<&[u8]> = haps.iter().map(|h| h.as_slice()).collect();
    let got = score_read_haps_logless(backend, &read, &quals, &hap_refs, &ins, &del, &gcp).unwrap();
    for (i, hap) in hap_refs.iter().enumerate() {
        let scalar = logless_pairhmm_likelihood(&read, &quals, hap, &ins, &del, &gcp).unwrap();
        assert!(
            pairhmm_fp_eq(got[i], scalar, policy) || (got[i].is_infinite() && scalar.is_infinite()),
            "f32 path i={i}: got={} scalar={}",
            got[i],
            scalar
        );
    }
}

#[test]
fn parse_pair_hmm_aliases() {
    use gatk_haplotypecaller::parse_pair_hmm_impl;
    assert_eq!(
        parse_pair_hmm_impl("scalar").unwrap(),
        PairHmmImpl::Log10PairHmm
    );
    assert_eq!(parse_pair_hmm_impl("avx").unwrap(), PairHmmImpl::Simd);
    assert_eq!(
        parse_pair_hmm_impl("FASTEST_AVAILABLE").unwrap(),
        PairHmmImpl::FastestAvailable
    );
}
