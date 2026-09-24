/// L13-B: PairHMM → biallelic GL / informative-AD score stage.
/// Owns marginalize + normalize + [`genotype_from_marginalized_rows`] (via the former
/// `genotype_from_allele_mapping` body). Called after [`SiteMap`] and before [`SiteReshape`].
pub(crate) struct SiteScore;

/// Typed ALT haplotype index pool for biallelic marginalize (L13-C4).
#[derive(Debug, Clone)]
pub(crate) struct AltHapSubset(pub Vec<HaplotypeIndex>);

impl AltHapSubset {
    pub(crate) fn as_slice(&self) -> &[HaplotypeIndex] {
        &self.0
    }
}

thread_local! {
    static LAST_ROWS_CACHE_HIT: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    static LAST_SITE_SCORE_INNER: std::cell::RefCell<Option<SiteScoreInnerTrace>> =
        const { std::cell::RefCell::new(None) };
    static LAST_ROWS_CACHE_LOOKUP: std::cell::RefCell<Option<RegionLikelihoodRowsLookupTrace>> =
        const { std::cell::RefCell::new(None) };
    static LAST_INSERT_SPARSE: std::cell::Cell<u64> = const { std::cell::Cell::new(0) };
    static LAST_INSERT_HAD: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    static ROWS_CACHE_LOOKUP_LOG: std::cell::RefCell<Vec<RegionLikelihoodRowsLookupTrace>> =
        const { std::cell::RefCell::new(Vec::new()) };
    /// 0 = production, 1 = disable (always rebuild), 2 = content-hash key.
    static ROWS_CACHE_DIAGNOSTIC: std::cell::Cell<u8> = const { std::cell::Cell::new(0) };
}

/// 6R.224 investigation-only: last `with_region_likelihood_rows` lookup key + input identity.
/// Production HIT/MISS rules are unchanged; this only records the last invocation.
#[derive(Debug, Clone, Copy)]
pub struct RegionLikelihoodRowsLookupTrace {
    pub hit: bool,
    pub input_ptr: usize,
    pub input_len: usize,
    pub n_haps: usize,
    pub input_sparse_hash: u64,
    pub input_read_set_hash: u64,
    pub input_n_unique_reads: usize,
    pub stored_ptr: Option<usize>,
    pub stored_len: Option<usize>,
    pub stored_n_haps: Option<usize>,
    pub stored_n_rows: Option<usize>,
    /// HIT and input sparse-cell hash equals the last inserted population.
    pub legitimate_hit: bool,
    /// HIT but input sparse-cell hash differs from the last inserted population.
    pub invalid_alias: bool,
}

fn hash_sparse_likelihoods(likelihoods: &[RegionReadLikelihood]) -> (u64, u64, usize) {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};
    let mut sparse = DefaultHasher::new();
    let mut reads = std::collections::BTreeSet::new();
    for rl in likelihoods {
        let ri = rl.read_index.get();
        ri.hash(&mut sparse);
        rl.haplotype_index.get().hash(&mut sparse);
        ((rl.log10_likelihood * 1e9).round() as i64).hash(&mut sparse);
        reads.insert(ri);
    }
    let mut rs = DefaultHasher::new();
    for ri in &reads {
        ri.hash(&mut rs);
    }
    (sparse.finish(), rs.finish(), reads.len())
}

pub(crate) fn record_region_likelihood_rows_lookup(
    key: (usize, usize, usize),
    hit: bool,
    likelihoods: &[RegionReadLikelihood],
    stored: Option<&(Vec<(usize, usize, u64)>, usize, Vec<ReadLikelihoodRow>)>,
) {
    let (sparse, read_set, n_unique) = hash_sparse_likelihoods(likelihoods);
    let last = LAST_INSERT_SPARSE.get();
    let had = LAST_INSERT_HAD.get();
    let legitimate_hit = hit && had && sparse == last;
    let invalid_alias = hit && had && sparse != last;
    if !hit {
        LAST_INSERT_SPARSE.set(sparse);
        LAST_INSERT_HAD.set(true);
    }
    let trace = RegionLikelihoodRowsLookupTrace {
        hit,
        input_ptr: key.0,
        input_len: key.1,
        n_haps: key.2,
        input_sparse_hash: sparse,
        input_read_set_hash: read_set,
        input_n_unique_reads: n_unique,
        stored_ptr: stored.map(|s| s.0.len()),
        stored_len: stored.map(|s| s.0.len()),
        stored_n_haps: stored.map(|s| s.1),
        stored_n_rows: stored.map(|s| s.2.len()),
        legitimate_hit,
        invalid_alias,
    };
    LAST_ROWS_CACHE_LOOKUP.with(|slot| {
        *slot.borrow_mut() = Some(trace);
    });
    ROWS_CACHE_LOOKUP_LOG.with(|log| {
        log.borrow_mut().push(trace);
    });
}

/// 6R.224 investigation-only: take the last [`with_region_likelihood_rows`] lookup snapshot.
#[doc(hidden)]
pub fn take_last_region_likelihood_rows_lookup_trace() -> Option<RegionLikelihoodRowsLookupTrace> {
    LAST_ROWS_CACHE_LOOKUP.with(|slot| slot.borrow_mut().take())
}

/// 6R.225 investigation-only: production=0, disable=1, content-hash key=2.
#[doc(hidden)]
pub fn set_region_likelihood_rows_cache_diagnostic(mode: u8) {
    ROWS_CACHE_DIAGNOSTIC.set(mode);
}

pub(crate) fn reset_rows_cache_insert_identity() {
    LAST_INSERT_HAD.set(false);
    LAST_INSERT_SPARSE.set(0);
    ROWS_CACHE_LOOKUP_LOG.with(|log| log.borrow_mut().clear());
}

/// 6R.225 investigation-only: drain lookup log (legitimate vs alias hits).
#[doc(hidden)]
pub fn take_region_likelihood_rows_lookup_log() -> Vec<RegionLikelihoodRowsLookupTrace> {
    ROWS_CACHE_LOOKUP_LOG.with(|log| std::mem::take(&mut *log.borrow_mut()))
}

/// Exact sparse-cell identity used as the production TLS cache key (6R.226).
/// Each cell is `(read_index, haplotype_index, log10_likelihood.to_bits())`.
/// This is equality of the production sparse population, not a hash.
pub(crate) fn sparse_population_identity(
    likelihoods: &[RegionReadLikelihood],
) -> Vec<(usize, usize, u64)> {
    likelihoods
        .iter()
        .map(|rl| {
            (
                rl.read_index.get(),
                rl.haplotype_index.get(),
                rl.log10_likelihood.to_bits(),
            )
        })
        .collect()
}

/// 6R.226 investigation-only: exact sparse-cell identity + `n_haps`.
#[doc(hidden)]
pub fn region_likelihood_rows_logical_identity(
    likelihoods: &[RegionReadLikelihood],
    n_haplotypes: usize,
) -> (Vec<(usize, usize, u64)>, usize) {
    (sparse_population_identity(likelihoods), n_haplotypes)
}

/// Diagnostic disable: rebuild without touching the production slot.
pub(crate) fn diagnostic_rows_cache_bypass(
    likelihoods: &[RegionReadLikelihood],
    n_haplotypes: usize,
) -> Option<Vec<ReadLikelihoodRow>> {
    if ROWS_CACHE_DIAGNOSTIC.get() != 1 {
        return None;
    }
    LAST_ROWS_CACHE_HIT.set(false);
    let rows = region_likelihoods_to_rows_uncached(likelihoods, n_haplotypes);
    record_region_likelihood_rows_lookup(
        (
            likelihoods.as_ptr() as usize,
            likelihoods.len(),
            n_haplotypes,
        ),
        false,
        likelihoods,
        None,
    );
    Some(rows)
}

pub(crate) fn diagnostic_rows_cache_key(
    likelihoods: &[RegionReadLikelihood],
    n_haplotypes: usize,
) -> (usize, usize, usize) {
    if ROWS_CACHE_DIAGNOSTIC.get() == 2 {
        let (sparse, _, _) = hash_sparse_likelihoods(likelihoods);
        return (sparse as usize, likelihoods.len(), n_haplotypes);
    }
    (
        likelihoods.as_ptr() as usize,
        likelihoods.len(),
        n_haplotypes,
    )
}

/// 6R.223 investigation-only: SiteScore inner snapshot (dense rows → 52×2 → AD/PL).
/// Production genotype math is unchanged; this only records the last invocation.
#[derive(Debug, Clone)]
pub struct SiteScoreInnerTrace {
    pub cache_hit: bool,
    pub dense_n_rows: usize,
    pub dense_n_cols: usize,
    pub dense_hash: u64,
    pub marg_n_rows: usize,
    pub marg_n_cols: usize,
    pub marg_hash: u64,
    /// Stable `(read_index, REF×1e9, ALT×1e9)` allele-likelihood cells.
    pub marg: Vec<(usize, i64, i64)>,
    pub ad: Vec<i32>,
    pub pl: Vec<i32>,
    /// Snapshot of the `with_region_likelihood_rows` lookup that fed this SiteScore.
    pub lookup: Option<RegionLikelihoodRowsLookupTrace>,
}

fn hash_ll_rows(rows: &[ReadLikelihoodRow]) -> u64 {
    use std::collections::hash_map::DefaultHasher;
    use std::hash::{Hash, Hasher};
    let mut h = DefaultHasher::new();
    for row in rows {
        row.read_index.hash(&mut h);
        for v in &row.haplotype_log10_likelihoods {
            ((v * 1e9).round() as i64).hash(&mut h);
        }
    }
    h.finish()
}

fn record_site_score_inner(
    rows: &[ReadLikelihoodRow],
    marg: &[ReadLikelihoodRow],
    gt: &RegionGenotypeResult,
) {
    let dense_cols = rows
        .first()
        .map(|r| r.haplotype_log10_likelihoods.len())
        .unwrap_or(0);
    let mut cells: Vec<(usize, i64, i64)> = marg
        .iter()
        .map(|row| {
            let lr = row
                .haplotype_log10_likelihoods
                .first()
                .copied()
                .unwrap_or(f64::NAN);
            let la = row
                .haplotype_log10_likelihoods
                .get(1)
                .copied()
                .unwrap_or(f64::NAN);
            (
                row.read_index,
                (lr * 1e9).round() as i64,
                (la * 1e9).round() as i64,
            )
        })
        .collect();
    cells.sort_unstable_by_key(|c| c.0);
    LAST_SITE_SCORE_INNER.with(|slot| {
        *slot.borrow_mut() = Some(SiteScoreInnerTrace {
            cache_hit: LAST_ROWS_CACHE_HIT.get(),
            dense_n_rows: rows.len(),
            dense_n_cols: dense_cols,
            dense_hash: hash_ll_rows(rows),
            marg_n_rows: marg.len(),
            marg_n_cols: marg
                .first()
                .map(|r| r.haplotype_log10_likelihoods.len())
                .unwrap_or(0),
            marg_hash: hash_ll_rows(marg),
            marg: cells,
            ad: gt.format.ad_as_i32(),
            pl: gt.format.pl_as_i32(),
            lookup: LAST_ROWS_CACHE_LOOKUP.with(|s| *s.borrow()),
        });
    });
}

/// 6R.223 investigation-only: take the last [`SiteScore::from_allele_mapping`] inner snapshot.
#[doc(hidden)]
pub fn take_last_site_score_inner_trace() -> Option<SiteScoreInnerTrace> {
    LAST_SITE_SCORE_INNER.with(|slot| slot.borrow_mut().take())
}

impl SiteScore {
    /// Score a site from allele↔haplotype mapping and region likelihood rows.
    pub(crate) fn from_allele_mapping(
        likelihoods: &[RegionReadLikelihood],
        haplotypes: &[Haplotype],
        mapping: &AlleleHaplotypeMapping,
        event: &VariationEvent,
        ref_bytes: &[u8],
        pad_start_1based: u64,
        max_mnp_distance: usize,
        contig: &str,
        config: &HcGenotypingConfig,
    ) -> GatkResult<RegionGenotypeResult> {
        let ref_hap = haplotypes
            .iter()
            .find(|h| h.is_reference)
            .or_else(|| haplotypes.first())
            .ok_or_else(|| {
                gatk_common::GatkError::algorithm(
                    "SiteScore::from_allele_mapping: haplotype list is empty",
                )
            })?;
        let profiling = crate::hc_profile::enabled();
        let t_marg = profiling.then(std::time::Instant::now);
        let out = with_region_likelihood_rows(likelihoods, haplotypes.len(), |rows| {
            let ref_pool = ref_hap_indices_for_genotype_marginalization(
                mapping,
                haplotypes,
                config,
                Some(event),
            );
            let alt_pool = AltHapSubset(alt_hap_indices_for_genotype_marginalization(
                mapping,
                haplotypes,
                event,
                ref_hap,
                pad_start_1based,
                ref_bytes,
                max_mnp_distance,
                contig,
                config,
            ));
            let mut marg =
                marginalize_rows_to_biallelic_alleles(rows, &ref_pool, alt_pool.as_slice());
            if config.enable_java_strict() {
                // Java AlleleLikelihoodMatrixMapper normalize applies at all strict sites (not only sparse).
                apply_java_marginal_normalize_gap(&mut marg);
            }
            if let Some(t0) = t_marg {
                crate::hc_profile::note_marginalize_wall(t0.elapsed());
            }
            let t_enum = profiling.then(std::time::Instant::now);
            let out = genotype_from_marginalized_rows(&marg, haplotypes, config);
            if let Some(t0) = t_enum {
                crate::hc_profile::note_genotype_enum_wall(t0.elapsed());
            }
            if let Ok(ref gt) = out {
                record_site_score_inner(rows, &marg, gt);
            }
            out
        });
        out
    }
}

/// Compatibility wrapper — production call sites may use [`SiteScore::from_allele_mapping`].
fn genotype_from_allele_mapping(
    likelihoods: &[RegionReadLikelihood],
    haplotypes: &[Haplotype],
    mapping: &AlleleHaplotypeMapping,
    event: &VariationEvent,
    ref_bytes: &[u8],
    pad_start_1based: u64,
    max_mnp_distance: usize,
    contig: &str,
    config: &HcGenotypingConfig,
) -> GatkResult<RegionGenotypeResult> {
    SiteScore::from_allele_mapping(
        likelihoods,
        haplotypes,
        mapping,
        event,
        ref_bytes,
        pad_start_1based,
        max_mnp_distance,
        contig,
        config,
    )
}
