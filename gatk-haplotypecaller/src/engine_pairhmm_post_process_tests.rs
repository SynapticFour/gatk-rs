use super::*;

#[test]
fn filter_drops_poorly_modeled_read_keeps_informative() {
    let mut ll = vec![
        RegionReadLikelihood {
            read_index: crate::bio_ids::ReadIndex::new(0),
            haplotype_index: crate::bio_ids::HaplotypeIndex::new(0),
            log10_likelihood: -50.0,
        },
        RegionReadLikelihood {
            read_index: crate::bio_ids::ReadIndex::new(0),
            haplotype_index: crate::bio_ids::HaplotypeIndex::new(1),
            log10_likelihood: -7.0,
        },
        RegionReadLikelihood {
            read_index: crate::bio_ids::ReadIndex::new(1),
            haplotype_index: crate::bio_ids::HaplotypeIndex::new(0),
            log10_likelihood: -49.0,
        },
        RegionReadLikelihood {
            read_index: crate::bio_ids::ReadIndex::new(1),
            haplotype_index: crate::bio_ids::HaplotypeIndex::new(1),
            log10_likelihood: -44.0,
        },
    ];
    normalize_region_read_likelihoods(&mut ll, &[0, 1]);
    let seq =
        b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    let qual =
        b"########################################################################################";
    let mut keep = rust_htslib::bam::Record::new();
    keep.set(b"ok", None, seq, qual);
    let mut drop = rust_htslib::bam::Record::new();
    drop.set(b"bad", None, seq, qual);
    let filtered = filter_poorly_modeled_region_read_likelihoods(
        &ll,
        &crate::shared_bam::share_records(vec![keep, drop]),
        None,
    );
    assert!(
        filtered.iter().any(|e| e.read_index.get() == 0),
        "informative read kept"
    );
    assert!(
        !filtered.iter().any(|e| e.read_index.get() == 1),
        "poorly modeled read dropped"
    );
}

#[test]
fn filter_drops_all_reads_when_none_pass_threshold() {
    let mut ll = vec![
        RegionReadLikelihood {
            read_index: crate::bio_ids::ReadIndex::new(0),
            haplotype_index: crate::bio_ids::HaplotypeIndex::new(0),
            log10_likelihood: -16.5526,
        },
        RegionReadLikelihood {
            read_index: crate::bio_ids::ReadIndex::new(0),
            haplotype_index: crate::bio_ids::HaplotypeIndex::new(1),
            log10_likelihood: -13.0042,
        },
    ];
    normalize_region_read_likelihoods(&mut ll, &[0, 1]);
    let mut rec = rust_htslib::bam::Record::new();
    rec.set(
        b"r2",
        None,
        b"ACGTACGTACGTACGTACGT",
        b"####################",
    );
    let reads = [crate::shared_bam::share_record(rec)];
    let filtered = filter_poorly_modeled_region_read_likelihoods(&ll, &reads, None);
    assert!(
        filtered.is_empty(),
        "P12 upstream marginal read class max LL -13.0 below Java static threshold -8.0"
    );
    let post = post_process_pairhmm_likelihoods(ll.clone(), &reads, &[], true, None);
    assert!(
        post.is_empty(),
        "Java does not retain full matrix when no read passes filterPoorlyModeledEvidence"
    );
}

/// Java `searchBestAllele` walks every likelihood-object haplotype column. A subset
/// mask (EventMap supporters) must not be used as `max_P` (6R.127/6R.128).
#[test]
fn normalize_max_p_uses_every_likelihood_column() {
    const CAP: f64 = -4.5;
    let mut subset = vec![
        RegionReadLikelihood {
            read_index: crate::bio_ids::ReadIndex::new(0),
            haplotype_index: crate::bio_ids::HaplotypeIndex::new(0),
            log10_likelihood: -10.0,
        },
        RegionReadLikelihood {
            read_index: crate::bio_ids::ReadIndex::new(0),
            haplotype_index: crate::bio_ids::HaplotypeIndex::new(1),
            log10_likelihood: -2.0,
        },
        RegionReadLikelihood {
            read_index: crate::bio_ids::ReadIndex::new(0),
            haplotype_index: crate::bio_ids::HaplotypeIndex::new(2),
            log10_likelihood: -20.0,
        },
    ];
    let mut all_cols = subset.clone();
    normalize_region_read_likelihoods(&mut subset, &[0]);
    normalize_region_read_likelihoods(&mut all_cols, &[0, 1, 2]);
    assert!(
        (subset[2].log10_likelihood - (-10.0 + CAP)).abs() < 1e-12,
        "subset mask floors against hap 0 only"
    );
    assert!(
        (all_cols[2].log10_likelihood - (-2.0 + CAP)).abs() < 1e-12,
        "full likelihood-object range floors against the true max column"
    );
    assert!(
        all_cols[1].log10_likelihood > subset[2].log10_likelihood,
        "including every column raises the floor vs an EventMap-style subset"
    );
}
