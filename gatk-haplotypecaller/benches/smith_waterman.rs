//! Microbenchmarks for Smith-Waterman on production-like haplotype/read lengths.
#![allow(clippy::result_large_err)]

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use gatk_haplotypecaller::smith_waterman::{align, SwOverhangStrategy, SwParameters};

use std::hint::black_box;
use std::time::Duration;

fn smoke_criterion() -> Criterion {
    Criterion::default()
        .warm_up_time(Duration::from_millis(300))
        .measurement_time(Duration::from_secs(2))
        .sample_size(40)
}

fn make_seqs(ref_len: usize, alt_len: usize, mismatch_every: usize) -> (Vec<u8>, Vec<u8>) {
    let bases = [b'A', b'C', b'G', b'T'];
    let reference: Vec<u8> = (0..ref_len).map(|i| bases[i % 4]).collect();
    let mut alternate: Vec<u8> = (0..alt_len).map(|i| bases[(i + 1) % 4]).collect();
    if mismatch_every > 0 {
        for i in (0..alt_len).step_by(mismatch_every) {
            alternate[i] = b'N';
        }
    }
    // Near-match interior so SoftClip still hits DP (not lastIndexOf) often.
    let copy = ref_len.min(alt_len).saturating_sub(8);
    if copy > 0 {
        alternate[..copy].copy_from_slice(&reference[..copy]);
    }
    (reference, alternate)
}

/// Pad like `haplotype_cigar::SW_PAD` (±10 N) around a haplotype-sized core.
fn make_padded_indel_pair(core_ref: usize, core_alt: usize) -> (Vec<u8>, Vec<u8>) {
    let pad = b"NNNNNNNNNN";
    let (r, a) = make_seqs(core_ref, core_alt, 11);
    let mut reference = Vec::with_capacity(pad.len() * 2 + r.len());
    reference.extend_from_slice(pad);
    reference.extend_from_slice(&r);
    reference.extend_from_slice(pad);
    let mut alternate = Vec::with_capacity(pad.len() * 2 + a.len());
    alternate.extend_from_slice(pad);
    alternate.extend_from_slice(&a);
    alternate.extend_from_slice(pad);
    (reference, alternate)
}

fn bench_sw_align(c: &mut Criterion) {
    let mut group = c.benchmark_group("smith_waterman_align");
    let params = SwParameters::gatk_haplotype_to_reference();
    // Legacy sizes + HC-typical: short hap, mid read×hap SoftClip, padded Indel.
    for &(ref_len, alt_len) in &[
        (64usize, 48usize),
        (128, 96),
        (256, 192),
        // Observed-ish: ~151-base read vs ~200-base haplotype (realign SoftClip).
        (200, 151),
        // Compact assembly haplotype pair before padding.
        (100, 100),
    ] {
        let id = format!("{ref_len}x{alt_len}");
        group.bench_with_input(BenchmarkId::new("soft_clip", &id), &id, |b, _| {
            let (reference, alternate) = make_seqs(ref_len, alt_len, 17);
            b.iter(|| {
                black_box(align(
                    black_box(&reference),
                    black_box(&alternate),
                    &params,
                    SwOverhangStrategy::SoftClip,
                ))
            });
        });
        group.bench_with_input(BenchmarkId::new("indel", &id), &id, |b, _| {
            let (reference, alternate) = make_seqs(ref_len, alt_len, 17);
            b.iter(|| {
                black_box(align(
                    black_box(&reference),
                    black_box(&alternate),
                    &params,
                    SwOverhangStrategy::Indel,
                ))
            });
        });
    }
    group.finish();
}

fn bench_sw_hc_workloads(c: &mut Criterion) {
    let mut group = c.benchmark_group("smith_waterman_hc");
    let hap_ref = SwParameters::gatk_haplotype_to_reference();
    let read_hap = SwParameters::gatk_read_to_best_haplotype();

    // Hap-to-ref Indel after ±10 N pad (assembly CIGAR path).
    for &core in &[80usize, 150, 250] {
        let (reference, alternate) = make_padded_indel_pair(core, core + 3);
        group.bench_with_input(BenchmarkId::new("padded_indel", core), &core, |b, _| {
            b.iter(|| {
                black_box(align(
                    black_box(&reference),
                    black_box(&alternate),
                    &hap_ref,
                    SwOverhangStrategy::Indel,
                ))
            });
        });
    }

    // Read→best-hap SoftClip (realign), mismatch so lastIndexOf misses.
    for &(hap_len, read_len) in &[(120usize, 100usize), (200, 151), (280, 151)] {
        let (hap, mut read) = make_seqs(hap_len, read_len, 0);
        // Force DP: break exact containment.
        if read_len > 10 {
            read[read_len / 2] = b'N';
        }
        let id = format!("hap{hap_len}_read{read_len}");
        group.bench_with_input(BenchmarkId::new("read_to_hap_soft", &id), &id, |b, _| {
            b.iter(|| {
                black_box(align(
                    black_box(&hap),
                    black_box(&read),
                    &read_hap,
                    SwOverhangStrategy::SoftClip,
                ))
            });
        });
    }

    // Exact substring fast path (common when read is contained in hap).
    {
        let hap: Vec<u8> = (0..200).map(|i| b"ACGT"[i % 4]).collect();
        let read = hap[20..170].to_vec();
        group.bench_function("exact_substring_fast_path", |b| {
            b.iter(|| {
                black_box(align(
                    black_box(&hap),
                    black_box(&read),
                    &read_hap,
                    SwOverhangStrategy::SoftClip,
                ))
            });
        });
    }

    group.finish();
}

criterion_group! {
    name = benches;
    config = smoke_criterion();
    targets = bench_sw_align, bench_sw_hc_workloads
}
criterion_main!(benches);
