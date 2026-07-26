//! Value-copy vs [`Arc<[u8]>`] sharing for padded active-region reference bytes.
//! Models the pre-refactor pattern (each pipeline stage `Vec::clone`s the window)
//! against the post-refactor pattern (`AssemblyResultSet::reference_bases_shared`
//! refcount bump only).
//! ```text
//! cargo bench -p gatk-haplotypecaller --bench shared_reference_arc
//! ```
#![allow(clippy::result_large_err)]

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use gatk_haplotypecaller::{
    AssemblyResultSet, AssemblyStatus, Cigar, CigarOperator, Haplotype, ThreadingAssemblyResult,
};
use std::hint::black_box;
use std::sync::Arc;

/// Medium synthetic padded region (~3 kb) — typical HC active-region window scale.
const REGION_BP: usize = 3072;
/// Stages that historically re-owned the window (trim / sync / supplement / genotype prep).
const PIPELINE_STAGES: usize = 12;

fn synthetic_ref(len: usize) -> Vec<u8> {
    let bases = [b'A', b'C', b'G', b'T'];
    (0..len).map(|i| bases[i % 4]).collect()
}

#[inline(never)]
fn touch_bytes(bytes: &[u8]) -> u64 {
    // Keep the optimizer from eliding copies (checksum depends on payload).
    let mut h = 0u64;
    for (i, &b) in bytes.iter().enumerate().step_by(17) {
        h = h
            .wrapping_mul(131)
            .wrapping_add(b as u64)
            .wrapping_add(i as u64);
    }
    h
}

/// Pre-audit pattern: each stage deep-copies the padded reference.
fn value_copy_pipeline(ref_bytes: &[u8], stages: usize) -> u64 {
    let mut owned = ref_bytes.to_vec();
    let mut acc = 0u64;
    for stage in 0..stages {
        let stage_copy = owned.clone();
        acc = acc
            .wrapping_add(touch_bytes(&stage_copy))
            .wrapping_add(stage as u64);
        owned = stage_copy;
    }
    acc
}

/// Post-audit pattern: stages share [`Arc<[u8]>`] (refcount only).
fn arc_share_pipeline(ref_bytes: Arc<[u8]>, stages: usize) -> u64 {
    let mut shared = ref_bytes;
    let mut acc = 0u64;
    for stage in 0..stages {
        let stage_share = Arc::clone(&shared);
        acc = acc
            .wrapping_add(touch_bytes(&stage_share))
            .wrapping_add(stage as u64);
        shared = stage_share;
    }
    acc
}

fn minimal_assembly_result_set(ref_bases: &[u8]) -> AssemblyResultSet {
    let mut ref_hap = Haplotype::new(ref_bases.to_vec(), true);
    let mut cigar = Cigar::new();
    cigar.push(ref_bases.len(), CigarOperator::Match);
    ref_hap.cigar = Some(cigar);
    let result = ThreadingAssemblyResult {
        status: AssemblyStatus::JustAssembledReference,
        haplotypes: vec![ref_hap],
        kmer_size: 10,
        event_maps: Vec::new(),
    };
    AssemblyResultSet::from_assembly_for_calling(&result, ref_bases, 1, "chr21", 0)
}

fn bench_pipeline_share(c: &mut Criterion) {
    let mut group = c.benchmark_group("shared_reference_arc");
    group.throughput(Throughput::Bytes((REGION_BP * PIPELINE_STAGES) as u64));

    let ref_vec = synthetic_ref(REGION_BP);
    let ref_arc: Arc<[u8]> = Arc::from(ref_vec.as_slice());

    group.bench_function(
        BenchmarkId::new(
            "value_copy_stages",
            format!("{REGION_BP}x{PIPELINE_STAGES}"),
        ),
        |b| b.iter(|| black_box(value_copy_pipeline(black_box(&ref_vec), PIPELINE_STAGES))),
    );
    group.bench_function(
        BenchmarkId::new("arc_share_stages", format!("{REGION_BP}x{PIPELINE_STAGES}")),
        |b| {
            b.iter(|| {
                black_box(arc_share_pipeline(
                    black_box(Arc::clone(&ref_arc)),
                    PIPELINE_STAGES,
                ))
            })
        },
    );

    // Fan-out only: models passing the window into many call sites that retain a handle
    // without re-scanning (the dominant cost avoided by Arc).
    group.bench_function(
        BenchmarkId::new("assembly_result_set_arc_fanout", REGION_BP),
        |b| {
            let result = minimal_assembly_result_set(&ref_vec);
            b.iter(|| {
                let mut handles = Vec::with_capacity(PIPELINE_STAGES);
                for _ in 0..PIPELINE_STAGES {
                    handles.push(result.reference_bases_shared());
                }
                black_box(handles)
            })
        },
    );
    group.bench_function(
        BenchmarkId::new("assembly_result_set_value_fanout", REGION_BP),
        |b| {
            let result = minimal_assembly_result_set(&ref_vec);
            let owned = result.reference_bases().to_vec();
            b.iter(|| {
                let mut handles = Vec::with_capacity(PIPELINE_STAGES);
                for _ in 0..PIPELINE_STAGES {
                    handles.push(owned.clone());
                }
                black_box(handles)
            })
        },
    );

    group.finish();
}

criterion_group!(benches, bench_pipeline_share);
criterion_main!(benches);
