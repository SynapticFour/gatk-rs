//! Microbenchmarks for Smith-Waterman (nested → flat DP matrices).
#![allow(clippy::result_large_err)]

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use gatk_haplotypecaller::smith_waterman::{align, SwOverhangStrategy, SwParameters};

fn make_seqs(ref_len: usize, alt_len: usize, mismatch_every: usize) -> (Vec<u8>, Vec<u8>) {
    let bases = [b'A', b'C', b'G', b'T'];
    let reference: Vec<u8> = (0..ref_len).map(|i| bases[i % 4]).collect();
    let mut alternate: Vec<u8> = (0..alt_len).map(|i| bases[(i + 1) % 4]).collect();
    if mismatch_every > 0 {
        for i in (0..alt_len).step_by(mismatch_every) {
            alternate[i] = b'N';
        }
    }
    // Keep most of alt as a near-match substring of ref so SoftClip path still hits DP often.
    let copy = ref_len.min(alt_len).saturating_sub(8);
    if copy > 0 {
        alternate[..copy].copy_from_slice(&reference[..copy]);
    }
    (reference, alternate)
}

fn bench_sw_align(c: &mut Criterion) {
    let mut group = c.benchmark_group("smith_waterman_align");
    let params = SwParameters::gatk_haplotype_to_reference();
    for &(ref_len, alt_len) in &[(64usize, 48usize), (128, 96), (256, 192)] {
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

criterion_group!(benches, bench_sw_align);
criterion_main!(benches);
