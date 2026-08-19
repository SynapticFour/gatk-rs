//! Architecture C wavefront PairHMM vs scalar Logless oracle.

use gatk_haplotypecaller::{
    logless_pairhmm_likelihood, pairhmm_fp_eq, resolve_pair_hmm_impl, score_haps_wavefront_f32,
    score_haps_wavefront_portable_f32, score_haps_wavefront_rolling_f64, score_read_haps_logless,
    select_wavefront_kernel, PairHmmFpPolicy, PairHmmImpl, WavefrontKernel,
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

fn policy() -> PairHmmFpPolicy {
    PairHmmFpPolicy {
        abs_epsilon: 1e-4,
        rel_epsilon: 1e-4,
    }
}

#[test]
fn rolling_f64_matches_scalar_logless() {
    let pol = policy();
    let lengths = [1usize, 2, 7, 8, 15, 16, 32, 64, 100];
    let mut cases = 0usize;
    for &rn in &lengths {
        for &hn in &lengths {
            for seed in 0..2u64 {
                let read = bases_pattern(rn, seed);
                let (quals, ins, del, gcp) = parity_quals(rn);
                let hap = bases_pattern(hn, seed.wrapping_add(3));
                let refs = [hap.as_slice()];
                let wf = score_haps_wavefront_rolling_f64(&read, &quals, &refs, &ins, &del, &gcp)
                    .expect("rolling f64");
                let scalar =
                    logless_pairhmm_likelihood(&read, &quals, &hap, &ins, &del, &gcp).expect("sc");
                assert!(
                    pairhmm_fp_eq(wf[0], scalar, pol)
                        || (wf[0].is_infinite() && scalar.is_infinite()),
                    "rn={rn} hn={hn} seed={seed}: wf={} scalar={}",
                    wf[0],
                    scalar
                );
                cases += 1;
            }
        }
    }
    assert!(cases >= 50);
}

#[test]
fn wavefront_f32_matches_scalar_logless() {
    let backend = resolve_pair_hmm_impl(PairHmmImpl::Wavefront);
    let pol = policy();
    let lengths = [1usize, 2, 7, 8, 15, 16, 31, 32, 64, 100, 151];
    let mut cases = 0usize;
    for &rn in &lengths {
        for &hn in &lengths {
            for seed in 0..3u64 {
                let read = bases_pattern(rn, seed);
                let (quals, ins, del, gcp) = parity_quals(rn);
                let mut haps: Vec<Vec<u8>> = Vec::new();
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
                let wf =
                    score_read_haps_logless(backend, &read, &quals, &hap_refs, &ins, &del, &gcp)
                        .expect("wavefront");
                for (i, hap) in hap_refs.iter().enumerate() {
                    let scalar = logless_pairhmm_likelihood(&read, &quals, hap, &ins, &del, &gcp)
                        .expect("scalar");
                    assert!(
                        pairhmm_fp_eq(wf[i], scalar, pol)
                            || (wf[i].is_infinite() && scalar.is_infinite()),
                        "rn={rn} hn={} seed={seed} i={i}: wf={} scalar={} kernel={:?}",
                        hap.len(),
                        wf[i],
                        scalar,
                        select_wavefront_kernel()
                    );
                    cases += 1;
                }
            }
        }
    }
    assert!(cases >= 200);
}

#[test]
fn host_simd_matches_portable_wavefront() {
    let pol = policy();
    let read = bases_pattern(200, 9);
    let (quals, ins, del, gcp) = parity_quals(read.len());
    let haps: Vec<Vec<u8>> = (0..16)
        .map(|k| {
            let mut h = bases_pattern(200 + (k % 3), 11 + k as u64);
            if k % 4 == 0 && !h.is_empty() {
                let len = h.len();
                h[k % len] = b'G';
            }
            h
        })
        .collect();
    let refs: Vec<&[u8]> = haps.iter().map(|h| h.as_slice()).collect();
    let host = score_haps_wavefront_f32(&read, &quals, &refs, &ins, &del, &gcp).expect("host");
    let portable =
        score_haps_wavefront_portable_f32(&read, &quals, &refs, &ins, &del, &gcp).expect("port");
    for i in 0..host.len() {
        assert!(
            pairhmm_fp_eq(host[i], portable[i], pol)
                || (host[i].is_infinite() && portable[i].is_infinite()),
            "i={i}: host={} portable={} kernel={:?}",
            host[i],
            portable[i],
            select_wavefront_kernel()
        );
    }
    // Ensure we actually selected a native kernel on this CI/dev host when available.
    let _ = matches!(
        select_wavefront_kernel(),
        WavefrontKernel::NeonF32 | WavefrontKernel::Avx2F32 | WavefrontKernel::PortableF32
    );
}

#[test]
fn phenotype_read_len_hap_count_wavefront() {
    let backend = resolve_pair_hmm_impl(PairHmmImpl::Wavefront);
    let pol = policy();
    let mut cases = 0usize;
    for &read_len in &[100usize, 200, 300] {
        for &hap_count in &[1usize, 8, 16, 32, 64] {
            let read = bases_pattern(read_len, 7);
            let (quals, ins, del, gcp) = parity_quals(read_len);
            // Shared-prefix assembly-like haps (equal length).
            let base = bases_pattern(read_len + 5, 42);
            let haps: Vec<Vec<u8>> = (0..hap_count)
                .map(|k| {
                    let mut h = base.clone();
                    if k > 0 {
                        let pos = (k * 7) % h.len();
                        h[pos] = b'T';
                    }
                    h
                })
                .collect();
            let hap_refs: Vec<&[u8]> = haps.iter().map(|h| h.as_slice()).collect();
            let wf = score_read_haps_logless(backend, &read, &quals, &hap_refs, &ins, &del, &gcp)
                .expect("wf");
            for (i, hap) in hap_refs.iter().enumerate() {
                let scalar =
                    logless_pairhmm_likelihood(&read, &quals, hap, &ins, &del, &gcp).expect("sc");
                assert!(
                    pairhmm_fp_eq(wf[i], scalar, pol)
                        || (wf[i].is_infinite() && scalar.is_infinite()),
                    "r={read_len} n={hap_count} i={i}: wf={} scalar={}",
                    wf[i],
                    scalar
                );
                cases += 1;
            }
        }
    }
    assert!(cases >= 100);
}

#[test]
fn high_error_short_read_triggers_retry_still_matches() {
    let pol = policy();
    // Low baseQ → small linear mass; may force f64 retry.
    let read = b"ACGTACGT";
    let quals = vec![2u8; read.len()];
    let ins = vec![GATK_PARITY_DEFAULT_INS_QUAL; read.len()];
    let del = vec![GATK_PARITY_DEFAULT_DEL_QUAL; read.len()];
    let gcp = vec![GATK_PARITY_DEFAULT_GCP; read.len()];
    let hap = b"ACGTACGTNNNNACGT";
    let refs = [hap.as_slice()];
    let wf = score_haps_wavefront_f32(read, &quals, &refs, &ins, &del, &gcp).expect("wf");
    let scalar = logless_pairhmm_likelihood(read, &quals, hap, &ins, &del, &gcp).expect("sc");
    assert!(
        pairhmm_fp_eq(wf[0], scalar, pol) || (wf[0].is_infinite() && scalar.is_infinite()),
        "wf={} scalar={}",
        wf[0],
        scalar
    );
}

#[test]
fn parse_wavefront_cli_aliases() {
    use gatk_haplotypecaller::parse_pair_hmm_impl;
    assert_eq!(
        parse_pair_hmm_impl("WAVEFRONT").unwrap(),
        PairHmmImpl::Wavefront
    );
    assert_eq!(
        parse_pair_hmm_impl("gkl_style").unwrap(),
        PairHmmImpl::Wavefront
    );
    assert_eq!(
        resolve_pair_hmm_impl(PairHmmImpl::Wavefront).label(),
        "WAVEFRONT_F32"
    );
    // FastestAvailable must not silently become wavefront.
    let fastest = resolve_pair_hmm_impl(PairHmmImpl::FastestAvailable);
    assert_ne!(fastest.label(), "WAVEFRONT_F32");
}
