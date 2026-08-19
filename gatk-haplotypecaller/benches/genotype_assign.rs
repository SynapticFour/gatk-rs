//! Adversarial genotype-assignment benches (dense multi-allelic + multi-pass AD).
//!
//! Reproduces the mega TRACE shape: many sites × several alts × deep coverage,
//! where PairHMM is paid once but genotyping rescans pileups per allele.
#![allow(clippy::result_large_err)]

use std::hint::black_box;
use std::time::Duration;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use gatk_haplotypecaller::bio_ids::HaplotypeIndex;
use gatk_haplotypecaller::event_map::VariationEvent;
use gatk_haplotypecaller::genome_loc::GenomePosition;
use gatk_haplotypecaller::genotyping::ReadLikelihoodRow;
use gatk_haplotypecaller::hc_genotyping_engine::{
    biallelic_genotype_log10_likelihoods_gatk, marginalize_rows_to_biallelic_alleles,
};
use gatk_haplotypecaller::read_event_discovery::{
    clear_ad_decode_cache, read_allele_depths_at_locus, read_allele_depths_at_locus_dedupe_qname,
};
use gatk_haplotypecaller::shared_bam::{share_record, SharedBamRecord};
use rust_htslib::bam::record::{Cigar as HtsCigar, CigarString, Record};

fn smoke_criterion() -> Criterion {
    Criterion::default()
        .warm_up_time(Duration::from_millis(300))
        .measurement_time(Duration::from_secs(2))
        .sample_size(30)
}

fn make_read(seq: &[u8], pos_1based: u64, qname: &[u8]) -> SharedBamRecord {
    let mut rec = Record::new();
    let cigar = CigarString::from(vec![HtsCigar::Match(seq.len() as u32)]);
    let qual: Vec<u8> = vec![30; seq.len()];
    rec.set(qname, Some(&cigar), seq, &qual);
    rec.set_pos(i64::try_from(pos_1based.saturating_sub(1)).unwrap_or(0));
    share_record(rec)
}

/// Dense coverage pileup: R reads overlapping one locus (mega-region shape).
fn dense_pileup(reads: usize, locus: u64) -> Vec<SharedBamRecord> {
    let mut out = Vec::with_capacity(reads);
    let bases = [b'A', b'C', b'G', b'T'];
    for i in 0..reads {
        let mut seq = vec![b'A'; 151];
        // Mix REF=A and ALT=C at the locus offset (pos = locus-50).
        let off = 50usize;
        seq[off] = if i % 3 == 0 { b'C' } else { bases[i % 4] };
        let qn = format!("r{i}").into_bytes();
        out.push(make_read(&seq, locus.saturating_sub(50), &qn));
    }
    out
}

fn snp_event(loc: u64, alt: &str) -> VariationEvent {
    VariationEvent {
        contig: "chr21".to_string(),
        start_1based: GenomePosition::new_1based(loc),
        end_1based: GenomePosition::new_1based(loc),
        ref_allele: "A".to_string(),
        alt_allele: alt.to_string(),
    }
}

/// Production-like multi-pass AD: same reads rescanned P times per allele × A alleles × S sites.
fn multipass_ad_region(reads: &[SharedBamRecord], sites: usize, alleles: usize, passes: usize) {
    clear_ad_decode_cache();
    let pad = 1u64;
    for s in 0..sites {
        let loc = 9_826_233 + s as u64;
        for a in 0..alleles {
            let alt = match a % 3 {
                0 => "C",
                1 => "G",
                _ => "T",
            };
            let event = snp_event(loc, alt);
            for p in 0..passes {
                if p % 2 == 0 {
                    black_box(read_allele_depths_at_locus(reads, &event, pad));
                } else {
                    black_box(read_allele_depths_at_locus_dedupe_qname(reads, &event, pad));
                }
            }
        }
    }
}

fn make_likelihood_rows(reads: usize, haps: usize) -> Vec<ReadLikelihoodRow> {
    (0..reads)
        .map(|i| ReadLikelihoodRow {
            read_index: i,
            read_id: format!("read_{i}"),
            haplotype_log10_likelihoods: (0..haps)
                .map(|h| -0.1 * (i as f64 + 1.0) * (h as f64 + 1.0))
                .collect(),
        })
        .collect()
}

fn bench_multipass_ad(c: &mut Criterion) {
    let mut group = c.benchmark_group("genotype_dense_ad");
    // Mega TRACE-ish: deep coverage, many sites, several alts, P≈6 passes (production).
    for &(reads, sites, alleles, passes) in &[
        (64usize, 20usize, 3usize, 6usize),
        (256, 50, 4, 6),
        (512, 100, 5, 8),
    ] {
        let pileup = dense_pileup(reads, 9_826_233);
        let id = format!("R{reads}_S{sites}_A{alleles}_P{passes}");
        group.bench_with_input(BenchmarkId::new("multipass_ad", &id), &id, |b, _| {
            b.iter(|| {
                multipass_ad_region(
                    black_box(&pileup),
                    black_box(sites),
                    black_box(alleles),
                    black_box(passes),
                )
            });
        });
    }
    group.finish();
}

fn bench_marginalize_and_pl(c: &mut Criterion) {
    let mut group = c.benchmark_group("genotype_marginalize_pl");
    for &(reads, haps, alleles_as_sites) in &[(128usize, 16usize, 8usize), (512, 32, 16)] {
        let rows = make_likelihood_rows(reads, haps);
        let ref_pool: Vec<HaplotypeIndex> = vec![HaplotypeIndex::new(0)];
        let id = format!("R{reads}_H{haps}_A{alleles_as_sites}");
        group.bench_with_input(BenchmarkId::new("marg_then_pl", &id), &id, |b, _| {
            b.iter(|| {
                for a in 0..alleles_as_sites {
                    let alt = HaplotypeIndex::new(1 + (a % (haps - 1)));
                    let marg =
                        marginalize_rows_to_biallelic_alleles(black_box(&rows), &ref_pool, &[alt]);
                    black_box(biallelic_genotype_log10_likelihoods_gatk(&marg, 0, 1));
                }
            });
        });
    }
    group.finish();
}

/// Multi-allelic reshape: TLS borrow cache vs uncached rebuild (A alleles).
fn bench_likelihood_reshape_reuse(c: &mut Criterion) {
    use gatk_haplotypecaller::bio_ids::{HaplotypeIndex, ReadIndex};
    use gatk_haplotypecaller::hc_genotyping_engine::{
        region_likelihoods_to_rows_uncached_pub, with_region_likelihood_rows,
    };
    use gatk_haplotypecaller::region_read_likelihood::RegionReadLikelihood;

    let mut group = c.benchmark_group("genotype_likelihood_reshape");
    for &(reads, haps, alleles) in &[(128usize, 16usize, 8usize), (512, 32, 16)] {
        let mut sparse = Vec::with_capacity(reads * haps);
        for r in 0..reads {
            for h in 0..haps {
                sparse.push(RegionReadLikelihood {
                    read_index: ReadIndex::new(r),
                    haplotype_index: HaplotypeIndex::new(h),
                    log10_likelihood: -0.1 * (r as f64 + 1.0) * (h as f64 + 1.0),
                });
            }
        }
        let id = format!("R{reads}_H{haps}_A{alleles}");
        group.bench_with_input(BenchmarkId::new("cached_borrow", &id), &id, |b, _| {
            b.iter(|| {
                for _ in 0..alleles {
                    black_box(with_region_likelihood_rows(
                        black_box(&sparse),
                        haps,
                        |rows| rows.len(),
                    ));
                }
            });
        });
        group.bench_with_input(BenchmarkId::new("uncached_rebuild", &id), &id, |b, _| {
            b.iter(|| {
                for _ in 0..alleles {
                    black_box(region_likelihoods_to_rows_uncached_pub(
                        black_box(&sparse),
                        haps,
                    ));
                }
            });
        });
    }
    group.finish();
}

criterion_group! {
    name = benches;
    config = smoke_criterion();
    targets = bench_multipass_ad, bench_marginalize_and_pl, bench_likelihood_reshape_reuse
}
criterion_main!(benches);
