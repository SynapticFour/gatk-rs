//! GATK-aligned activity profile scaffolding (Phase 4, step 51).
//! Mirrors the data model and numeric thresholds from GATK’s
//! `org.broadinstitute.hellbender.utils.activityprofile` package and
//! `AssemblyRegionArgumentCollection` defaults used by
//! [`BandPassActivityProfile`](https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/utils/activityprofile/BandPassActivityProfile.java).

use gatk_common::{AssemblyRegionConfig, GatkError};
use std::sync::Arc;

/// GATK `AssemblyRegionArgumentCollection.DEFAULT_ACTIVE_PROB_THRESHOLD` (`--active-probability-threshold`).
pub const GATK_DEFAULT_ACTIVE_PROB_THRESHOLD: f64 =
    AssemblyRegionConfig::GATK_DEFAULT_ACTIVE_PROB_THRESHOLD;

/// GATK `AssemblyRegionArgumentCollection.DEFAULT_MAX_PROB_PROPAGATION_DISTANCE` (`--max-prob-propagation-distance`).
pub const GATK_DEFAULT_MAX_PROB_PROPAGATION_DISTANCE: u32 =
    AssemblyRegionConfig::GATK_DEFAULT_MAX_PROB_PROPAGATION_DISTANCE;

/// GATK `BandPassActivityProfile.MAX_FILTER_SIZE`.
pub const GATK_BAND_PASS_MAX_FILTER_SIZE: u32 =
    AssemblyRegionConfig::GATK_ACTIVITY_PROFILE_MAX_FILTER_SIZE;

/// GATK `BandPassActivityProfile.DEFAULT_SIGMA` (Gaussian kernel σ for the band-pass smoother).
pub const GATK_BAND_PASS_DEFAULT_SIGMA: f64 = AssemblyRegionConfig::GATK_ACTIVITY_PROFILE_SIGMA;

/// GATK `BandPassActivityProfile.MIN_PROB_TO_KEEP_IN_FILTER`.
pub const GATK_BAND_PASS_MIN_PROB_TO_KEEP_IN_FILTER: f64 = 1e-5;

#[inline]
fn root_two_pi() -> f64 {
    std::f64::consts::TAU.sqrt()
}

/// Evidence carried by one activity-profile locus (GATK `ActivityProfileState.Type`).
/// # Compiler-enforced
/// Soft-clip length exists only on [`ActivityEvidence::HighQualitySoftClips`] — impossible to pair
/// a plain locus with a clip length or an HQ-soft-clip locus without one.
/// # Java equivalence
/// GATK `ActivityProfileState.Type` (`NONE`, `HIGH_QUALITY_SOFT_CLIPS`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ActivityEvidence {
    Plain,
    HighQualitySoftClips { clip_bases: u32 },
}

/// Back-compat alias for dump / older call sites.
pub type ActivityProfileStateKind = ActivityEvidence;

/// One locus in the per-base activity profile (`ActivityProfileState` in GATK).
/// # Invariants
/// `pos` is a **1-based** single-base reference locus on `contig`.
/// Active coloring uses `active_prob > active_prob_threshold` ([`Self::is_active`]).
/// # Ownership
/// Contig is [`Arc<str>`] so band-pass / HQ expansion can clone loci without reallocating the name.
/// # Mutation
/// Band-pass profiles merge/replace states in a dense vector; individual states are typically replaced, not field-mutated.
/// # Biological assumptions
/// `active_prob` is probability mass that the locus warrants local assembly.
/// # Java equivalence
/// GATK `ActivityProfileState` (`org.broadinstitute.hellbender.utils.activityprofile`).
#[derive(Debug, Clone)]
pub struct ActivityProfileState {
    pub contig: Arc<str>,
    /// 1-based reference position (single-base locus).
    pub pos: u64,
    /// Smoothed / propagated probability mass for “active” at this locus.
    pub active_prob: f64,
    /// Pre-smoothing score (used by some evaluators, e.g. DRAGEN-GATK heuristics).
    pub original_active_prob: f64,
    pub evidence: ActivityEvidence,
}

impl ActivityProfileState {
    /// Plain active state (`Type.NONE`).
    /// Prefer passing an existing [`Arc<str>`] (cheap clone) over `&str` in hot loops.
    pub fn new(contig: impl Into<Arc<str>>, pos: u64, active_prob: f64) -> Self {
        Self {
            contig: contig.into(),
            pos,
            active_prob,
            original_active_prob: 0.0,
            evidence: ActivityEvidence::Plain,
        }
    }

    /// Soft-clip–derived state; `clip_bases` is capped later by max propagation distance (GATK `processState`).
    pub fn high_quality_soft_clips(
        contig: impl Into<Arc<str>>,
        pos: u64,
        active_prob: f64,
        clip_bases: u32,
    ) -> Self {
        Self {
            contig: contig.into(),
            pos,
            active_prob,
            original_active_prob: 0.0,
            evidence: ActivityEvidence::HighQualitySoftClips { clip_bases },
        }
    }

    /// GATK `Type` discriminant for dumps / legacy comparisons.
    #[inline]
    pub fn kind(&self) -> ActivityEvidence {
        self.evidence
    }

    /// HQ soft-clip length when present (GATK `resultValue`).
    #[inline]
    pub fn hq_soft_clip_bases(&self) -> Option<u32> {
        match self.evidence {
            ActivityEvidence::HighQualitySoftClips { clip_bases } => Some(clip_bases),
            ActivityEvidence::Plain => None,
        }
    }

    /// Offset in bp from a region start locus (GATK `ActivityProfileState#getOffset`).
    pub fn offset_from_region_start(&self, region_start_1based: u64) -> i64 {
        self.pos as i64 - region_start_1based as i64
    }

    /// GATK colors a locus “active” when `activeProb > activeProbThreshold`.
    pub fn is_active(&self, active_prob_threshold: f64) -> bool {
        self.active_prob > active_prob_threshold
    }
}

/// Band-pass construction parameters matching `AssemblyRegionIterator` → `BandPassActivityProfile` wiring.
/// # Invariants
/// `sigma > 0` for Gaussian kernel construction; adaptive filter trims below `MIN_PROB_TO_KEEP_IN_FILTER`.
/// Effective propagation distance = user max + resolved filter size.
/// # Ownership
/// Cloneable config snapshot; no internal allocation until kernel is built.
/// # Mutation
/// Immutable during a profile lifetime; callers replace the whole params bundle.
/// # Biological assumptions
/// Smooths sparse locus activity into contiguous active/inactive regions for assembly.
/// # Java equivalence
/// GATK `AssemblyRegionArgumentCollection` + `BandPassActivityProfile` filter/σ defaults.
/// Strictly positive σ for the band-pass Gaussian kernel.
/// # Compiler-enforced
/// [`PositiveSigma::try_new`] rejects `≤ 0`; only positive values construct.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct PositiveSigma(f64);

impl PositiveSigma {
    /// GATK `BandPassActivityProfile.DEFAULT_SIGMA` (`17.0`).
    pub const GATK_DEFAULT: Self = Self(GATK_BAND_PASS_DEFAULT_SIGMA);

    #[inline]
    pub fn try_new(sigma: f64) -> Option<Self> {
        if sigma > 0.0 && sigma.is_finite() {
            Some(Self(sigma))
        } else {
            None
        }
    }

    #[inline]
    pub const fn get(self) -> f64 {
        self.0
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct BandPassActivityProfileParams {
    pub max_prob_propagation_distance: u32,
    pub active_prob_threshold: f64,
    pub max_filter_size: u32,
    /// Gaussian kernel σ (`> 0`); type-enforced.
    pub sigma: PositiveSigma,
    pub adaptive_filter_size: bool,
}

impl BandPassActivityProfileParams {
    /// Same defaults as GATK `AssemblyRegionIterator` (adaptive kernel size, σ = 17, cap 50).
    pub fn gatk_haplotype_caller_defaults() -> Self {
        let cfg = AssemblyRegionConfig::default();
        Self {
            max_prob_propagation_distance: cfg.max_prob_propagation_distance,
            active_prob_threshold: cfg.active_prob_threshold,
            max_filter_size: cfg.activity_profile_max_filter_size,
            sigma: PositiveSigma::GATK_DEFAULT,
            adaptive_filter_size: true,
        }
    }

    /// Build from [`AssemblyRegionConfig`]; rejects non-positive / non-finite σ.
    pub fn try_from_assembly_region(cfg: &AssemblyRegionConfig) -> Result<Self, GatkError> {
        let sigma = PositiveSigma::try_new(cfg.activity_profile_sigma).ok_or_else(|| {
            GatkError::invalid_configuration(
                "activity_profile_sigma",
                format!(
                    "band-pass sigma must be finite and > 0, got {}",
                    cfg.activity_profile_sigma
                ),
            )
        })?;
        Ok(Self {
            max_prob_propagation_distance: cfg.max_prob_propagation_distance,
            active_prob_threshold: cfg.active_prob_threshold,
            max_filter_size: cfg.activity_profile_max_filter_size,
            sigma,
            adaptive_filter_size: true,
        })
    }

    /// Fallible constructor; rejects non-positive / non-finite σ.
    pub fn try_new(
        max_prob_propagation_distance: u32,
        active_prob_threshold: f64,
        max_filter_size: u32,
        sigma: f64,
        adaptive_filter_size: bool,
    ) -> Result<Self, GatkError> {
        let sigma = PositiveSigma::try_new(sigma).ok_or_else(|| {
            GatkError::invalid_argument(
                "sigma",
                format!("band-pass sigma must be finite and > 0, got {sigma}"),
            )
        })?;
        Ok(Self {
            max_prob_propagation_distance,
            active_prob_threshold,
            max_filter_size,
            sigma,
            adaptive_filter_size,
        })
    }

    /// Effective propagation distance for a band-pass profile (`BandPassActivityProfile#getMaxProbPropagationDistance`).
    pub fn effective_max_prob_propagation_distance(&self, resolved_filter_size: u32) -> u32 {
        self.max_prob_propagation_distance
            .saturating_add(resolved_filter_size)
    }

    /// Resolved half-width (Java `filterSize`) after optional adaptive trimming.
    pub fn resolved_filter_size(&self) -> u32 {
        self.resolved_filter_and_kernel().0
    }

    /// Normalized Gaussian kernel of length `2 * resolved_filter_size + 1` (Java `makeKernel` + `normalizeSumToOne`).
    pub fn normalized_kernel(&self) -> Vec<f64> {
        self.resolved_filter_and_kernel().1
    }

    /// Compute adaptive filter size and kernel together (one full-kernel build when adaptive).
    pub fn resolved_filter_and_kernel(&self) -> (u32, Vec<f64>) {
        if self.adaptive_filter_size {
            let full = make_gaussian_kernel(self.max_filter_size, self.sigma.get());
            let fs = adaptive_filter_size(&full, GATK_BAND_PASS_MIN_PROB_TO_KEEP_IN_FILTER);
            if fs == self.max_filter_size {
                (fs, full)
            } else {
                (fs, make_gaussian_kernel(fs, self.sigma.get()))
            }
        } else {
            let fs = self.max_filter_size;
            (fs, make_gaussian_kernel(fs, self.sigma.get()))
        }
    }
}

/// Unnormalized Gaussian PDF at `x` with mean `mean` and standard deviation `sd` (GATK `MathUtils.normalDistribution`).
pub fn normal_distribution(mean: f64, sd: f64, x: f64) -> f64 {
    debug_assert!(sd > 0.0, "sigma must be > 0 for a proper Gaussian kernel");
    (-((x - mean) * (x - mean)) / (2.0 * sd * sd)).exp() / (sd * root_two_pi())
}

/// GATK `BandPassActivityProfile.makeKernel` + `MathUtils.normalizeSumToOne`.
pub fn make_gaussian_kernel(filter_size: u32, sigma: f64) -> Vec<f64> {
    assert!(sigma > 0.0, "sigma must be > 0");
    let fs = filter_size as usize;
    let band = 2 * fs + 1;
    let mean = filter_size as f64;
    let mut kernel: Vec<f64> = (0..band)
        .map(|iii| normal_distribution(mean, sigma, iii as f64))
        .collect();
    let sum: f64 = kernel.iter().sum();
    debug_assert!(sum > 0.0, "kernel sum must be positive");
    for k in &mut kernel {
        *k /= sum;
    }
    kernel
}

/// GATK `BandPassActivityProfile.determineFilterSize`.
pub fn adaptive_filter_size(kernel: &[f64], min_prob_to_keep: f64) -> u32 {
    if kernel.is_empty() {
        return 0;
    }
    let middle = (kernel.len() - 1) / 2;
    let mut filter_end = middle;
    while filter_end > 0 && kernel[filter_end - 1] >= min_prob_to_keep {
        filter_end -= 1;
    }
    (middle - filter_end) as u32
}

/// GATK `ActivityProfile#processState` (HQ soft-clip expansion → plain loci; otherwise singleton).
pub fn activity_profile_base_process_state(
    just: &ActivityProfileState,
    max_prob_propagation_distance_with_filter: u32,
    contig_len: u64,
) -> Vec<ActivityProfileState> {
    match just.evidence {
        ActivityEvidence::HighQualitySoftClips { clip_bases } => {
            let num = clip_bases.min(max_prob_propagation_distance_with_filter);
            let mut out = Vec::new();
            for di in -(num as i64)..=(num as i64) {
                let p = (just.pos as i64).saturating_add(di);
                if p < 1 || p as u64 > contig_len {
                    continue;
                }
                // CLONE: needed — cheap `Arc<str>` bump; contig shared across expanded loci.
                out.push(ActivityProfileState::new(
                    just.contig.clone(),
                    p as u64,
                    just.active_prob,
                ));
            }
            out
        }
        ActivityEvidence::Plain => vec![ActivityProfileState::new(
            just.contig.clone(), // CLONE: needed — cheap `Arc<str>` bump
            just.pos,
            just.active_prob,
        )],
    }
}

/// GATK `BandPassActivityProfile#processState` (Gaussian window around the **input** locus).
pub fn band_pass_process_state(
    just_added: &ActivityProfileState,
    super_states: &[ActivityProfileState],
    filter_size: u32,
    kernel: &[f64],
    contig_len: u64,
) -> Vec<ActivityProfileState> {
    let fs = filter_size as usize;
    debug_assert_eq!(kernel.len(), 2 * fs + 1);
    let mut out = Vec::new();
    for super_s in super_states {
        if super_s.active_prob > 0.0 {
            for i in -(fs as i64)..=(fs as i64) {
                let start = (just_added.pos as i64).saturating_add(i);
                if start < 1 || start as u64 > contig_len {
                    continue;
                }
                let w = kernel[(i + fs as i64) as usize];
                let p = super_s.active_prob * w;
                // CLONE: needed — cheap `Arc<str>` bump for band-pass window states.
                out.push(ActivityProfileState::new(
                    just_added.contig.clone(),
                    start as u64,
                    p,
                ));
            }
        } else {
            // CLONE: needed — state is owned by the output vector (Arc contig + Copy fields).
            out.push(just_added.clone());
        }
    }
    out
}

/// Band-pass smoothed per-locus activity (`BandPassActivityProfile` + `ActivityProfile#add` merge semantics).
/// Input loci must arrive in **strictly increasing 1-based** order (GATK pileup order). Smoothed mass at
/// coordinates left of the first locus is dropped (same as GATK `incorporateSingleState` ignoring `position < 0`).
/// # Invariants
/// Dense `states[k]` sits at `region_start + k` on `contig`.
/// Input loci must be strictly increasing; gaps are filled per GATK merge rules.
/// # Ownership
/// Contig is [`Arc<str>`] shared into derived states; owns kernel and dense state vector.
/// # Mutation
/// `add` / process-state paths mutate `states` and tracking fields in place; region pops drain ready spans.
/// # Biological assumptions
/// Activity mass reflects variant evidence before cutting into assembly regions.
/// # Java equivalence
/// GATK `BandPassActivityProfile` / `ActivityProfile#add`.
#[derive(Debug, Clone)]
pub struct BandPassActivityProfile {
    contig: Arc<str>,
    contig_len: u64,
    params: BandPassActivityProfileParams,
    kernel: Vec<f64>,
    filter_size: usize,
    /// First locus added (profile coordinate origin).
    region_start: Option<u64>,
    /// Last **input** locus (pre-smoothing) for contiguity checks.
    last_input_pos: Option<u64>,
    /// Dense profile: `states[k]` sits at `region_start + k` on `contig`.
    states: Vec<ActivityProfileState>,
}

/// Region emitted from `ActivityProfile` cutting/merging.
/// # Invariants
/// Unpadded `start`/`end` and padded bounds are **1-based inclusive**, clipped to contig.
/// `padded_*` expands by `extension` (assembly-region padding).
/// # Ownership
/// Owns contig name; converted into [`crate::assembly_region_iterator::AssemblyRegion`] for apply.
/// # Mutation
/// Immutable cut result; padding applied at emission time.
/// # Biological assumptions
/// Active regions warrant local assembly; inactive may still get reference-confidence modeling.
/// # Java equivalence
/// GATK activity-profile cut regions feeding `AssemblyRegion` construction.
#[derive(Debug, Clone, PartialEq)]
pub struct ActivityProfileRegion {
    pub contig: String,
    /// Unpadded active/inactive region start, 1-based inclusive.
    pub start: u64,
    /// Unpadded active/inactive region end, 1-based inclusive.
    pub end: u64,
    pub is_active: bool,
    /// Padded bounds (`assemblyRegionExtension`), clipped to contig bounds.
    pub padded_start: u64,
    pub padded_end: u64,
    pub extension: u32,
}

impl BandPassActivityProfile {
    pub fn new(
        contig: impl Into<Arc<str>>,
        contig_len: u64,
        params: BandPassActivityProfileParams,
    ) -> Self {
        let (fs, kernel) = params.resolved_filter_and_kernel();
        let fs = fs as usize;
        debug_assert_eq!(kernel.len(), 2 * fs + 1);
        Self {
            contig: contig.into(),
            contig_len,
            params,
            kernel,
            filter_size: fs,
            region_start: None,
            last_input_pos: None,
            states: Vec::new(),
        }
    }

    #[inline]
    fn max_prob_propagation_distance_with_filter(&self) -> u32 {
        self.params
            .max_prob_propagation_distance
            .saturating_add(self.filter_size as u32)
    }

    /// Contiguous next locus (GATK `ActivityProfile#add` ordering).
    pub fn add(&mut self, state: ActivityProfileState) -> Result<(), GatkError> {
        if state.contig != self.contig {
            return Err(GatkError::argument(format!(
                "Activity profile contig mismatch: expected {}, got {}",
                self.contig, state.contig
            )));
        }
        if self.region_start.is_none() {
            self.region_start = Some(state.pos);
            self.last_input_pos = Some(state.pos);
            let derived = self.expand_smoothed(&state);
            self.incorporate_sorted(derived)?;
            return Ok(());
        }
        let last = self.last_input_pos.ok_or_else(|| {
            GatkError::argument("BandPassActivityProfile: missing last_input_pos")
        })?;
        if state.pos != last + 1 {
            return Err(GatkError::argument(format!(
                "Bad add to activity profile: locus {} not immediately after last {}",
                state.pos, last
            )));
        }
        self.last_input_pos = Some(state.pos);
        let derived = self.expand_smoothed(&state);
        self.incorporate_sorted(derived)?;
        Ok(())
    }

    fn expand_smoothed(&self, just: &ActivityProfileState) -> Vec<ActivityProfileState> {
        let max_hq = self.max_prob_propagation_distance_with_filter();
        let supers = activity_profile_base_process_state(just, max_hq, self.contig_len);
        band_pass_process_state(
            just,
            &supers,
            self.filter_size as u32,
            &self.kernel,
            self.contig_len,
        )
    }

    /// Merge derived loci into the dense profile (GATK `ActivityProfile#incorporateSingleState`).
    fn incorporate_sorted(
        &mut self,
        mut derived: Vec<ActivityProfileState>,
    ) -> Result<(), GatkError> {
        derived.sort_by_key(|s| s.pos);
        let mut merged: Vec<ActivityProfileState> = Vec::new();
        for s in derived {
            if let Some(t) = merged.last_mut() {
                if t.pos == s.pos {
                    t.active_prob += s.active_prob;
                    continue;
                }
            }
            merged.push(s);
        }

        let rs = self
            .region_start
            .expect("region_start set before incorporate_sorted");
        for s in merged {
            let offset = (s.pos as i64) - (rs as i64);
            if offset < 0 {
                continue;
            }
            let offset = offset as usize;
            if offset < self.states.len() {
                self.states[offset].active_prob += s.active_prob;
            } else if offset == self.states.len() {
                self.states.push(ActivityProfileState::new(
                    self.contig.clone(), // CLONE: needed — cheap `Arc<str>` bump into new state
                    s.pos,
                    s.active_prob,
                ));
            } else {
                return Err(GatkError::argument(format!(
                    "Activity profile gap: cannot add locus {} at offset {} (len {})",
                    s.pos,
                    offset,
                    self.states.len()
                )));
            }
        }
        Ok(())
    }

    pub fn contig(&self) -> &str {
        &self.contig
    }

    pub fn region_start(&self) -> Option<u64> {
        self.region_start
    }

    pub fn last_input_pos(&self) -> Option<u64> {
        self.last_input_pos
    }

    pub fn len(&self) -> usize {
        self.states.len()
    }

    pub fn is_empty(&self) -> bool {
        self.states.is_empty()
    }

    pub fn states(&self) -> &[ActivityProfileState] {
        &self.states
    }

    /// Sum of stored smoothed probabilities (useful invariant: one-hot input sums ≈1 per locus before merging).
    pub fn total_smoothed_mass(&self) -> f64 {
        self.states.iter().map(|s| s.active_prob).sum()
    }

    /// GATK `ActivityProfile#popNextReadyAssemblyRegion` — try one pop; call after each fed locus with
    /// `force_conversion = false`, then with `true` at interval boundaries (see `AssemblyRegionIterator`).
    pub fn try_pop_next_ready_region(
        &mut self,
        assembly_region_extension: u32,
        min_region_size: u32,
        max_region_size: u32,
        force_conversion: bool,
    ) -> Result<Option<ActivityProfileRegion>, GatkError> {
        if min_region_size == 0 || max_region_size == 0 {
            return Err(GatkError::argument(
                "min/max region size must be >= 1 for try_pop_next_ready_region",
            ));
        }
        self.pop_next_ready_region(
            assembly_region_extension,
            min_region_size,
            max_region_size,
            force_conversion,
        )
    }

    /// GATK `ActivityProfile#popReadyAssemblyRegions`.
    pub fn pop_ready_regions(
        &mut self,
        assembly_region_extension: u32,
        min_region_size: u32,
        max_region_size: u32,
        force_conversion: bool,
    ) -> Result<Vec<ActivityProfileRegion>, GatkError> {
        if min_region_size == 0 || max_region_size == 0 {
            return Err(GatkError::argument(
                "min/max region size must be >= 1 for pop_ready_regions",
            ));
        }
        let mut out = Vec::new();
        while let Some(r) = self.try_pop_next_ready_region(
            assembly_region_extension,
            min_region_size,
            max_region_size,
            force_conversion,
        )? {
            out.push(r);
        }
        Ok(out)
    }

    fn pop_next_ready_region(
        &mut self,
        assembly_region_extension: u32,
        min_region_size: u32,
        max_region_size: u32,
        force_conversion: bool,
    ) -> Result<Option<ActivityProfileRegion>, GatkError> {
        if self.states.is_empty() {
            return Ok(None);
        }

        // GATK force-flush behavior: trim states that lie beyond the current contiguous input span.
        if force_conversion {
            if let (Some(rs), Some(last_input)) = (self.region_start, self.last_input_pos) {
                let keep = last_input.saturating_sub(rs).saturating_add(1) as usize;
                if keep < self.states.len() {
                    self.states.truncate(keep);
                }
            }
        }
        if self.states.is_empty() {
            self.region_start = None;
            self.last_input_pos = None;
            return Ok(None);
        }

        let first = &self.states[0];
        let is_active_region = first.active_prob > self.params.active_prob_threshold;
        let maybe_end = self.find_end_of_region(
            is_active_region,
            min_region_size as usize,
            max_region_size as usize,
            force_conversion,
        );
        let end_idx = match maybe_end {
            Some(v) => v,
            None => return Ok(None),
        };

        let region_start = first.pos;
        let region_end = first.pos + end_idx as u64;
        let padded_start = region_start
            .saturating_sub(assembly_region_extension as u64)
            .max(1);
        let padded_end = region_end
            .saturating_add(assembly_region_extension as u64)
            .min(self.contig_len);

        self.states.drain(0..=end_idx);
        if self.states.is_empty() {
            self.region_start = None;
            self.last_input_pos = None;
        } else {
            self.region_start = Some(self.states[0].pos);
        }

        Ok(Some(ActivityProfileRegion {
            contig: String::from(self.contig.as_ref()),
            start: region_start,
            end: region_end,
            is_active: is_active_region,
            padded_start,
            padded_end,
            extension: assembly_region_extension,
        }))
    }

    fn find_end_of_region(
        &self,
        is_active_region: bool,
        min_region_size: usize,
        max_region_size: usize,
        force_conversion: bool,
    ) -> Option<usize> {
        if !force_conversion
            && self.states.len()
                < max_region_size + (self.max_prob_propagation_distance_with_filter() as usize)
        {
            return None;
        }

        let mut end_exclusive =
            self.find_first_activity_boundary(is_active_region, max_region_size);
        if is_active_region && end_exclusive == max_region_size {
            end_exclusive = self.find_best_cut_site(end_exclusive, min_region_size);
        }
        Some(end_exclusive - 1)
    }

    fn find_best_cut_site(&self, end_exclusive: usize, min_region_size: usize) -> usize {
        let mut min_i = end_exclusive - 1;
        let mut min_p = f64::MAX;
        let min_inclusive = min_region_size.saturating_sub(1);
        for i in (min_inclusive..=min_i).rev() {
            let cur = self.states[i].active_prob;
            if cur < min_p && self.is_minimum(i) {
                min_p = cur;
                min_i = i;
            }
        }
        min_i + 1
    }

    fn find_first_activity_boundary(
        &self,
        is_active_region: bool,
        max_region_size: usize,
    ) -> usize {
        let n_states = self.states.len();
        let mut end = 0usize;
        while end < n_states && end < max_region_size {
            if (self.states[end].active_prob > self.params.active_prob_threshold)
                != is_active_region
            {
                break;
            }
            end += 1;
        }
        end
    }

    fn is_minimum(&self, index: usize) -> bool {
        if index == self.states.len() - 1 || index < 1 {
            return false;
        }
        let p = self.states[index].active_prob;
        p <= self.states[index + 1].active_prob && p < self.states[index - 1].active_prob
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gatk_threshold_constants_match_assembly_region_defaults() {
        let cfg = AssemblyRegionConfig::default();
        assert_eq!(
            cfg.active_prob_threshold,
            GATK_DEFAULT_ACTIVE_PROB_THRESHOLD
        );
        assert_eq!(
            cfg.max_prob_propagation_distance,
            GATK_DEFAULT_MAX_PROB_PROPAGATION_DISTANCE
        );
        assert_eq!(
            cfg.activity_profile_max_filter_size,
            GATK_BAND_PASS_MAX_FILTER_SIZE
        );
        assert_eq!(cfg.activity_profile_sigma, GATK_BAND_PASS_DEFAULT_SIGMA);
    }

    #[test]
    fn band_pass_defaults_match_assembly_region_iterator() {
        let p = BandPassActivityProfileParams::gatk_haplotype_caller_defaults();
        assert!(p.adaptive_filter_size);
        assert_eq!(p.resolved_filter_size(), 50);
        let k = p.normalized_kernel();
        assert_eq!(k.len(), 101);
        let s: f64 = k.iter().sum();
        assert!((s - 1.0).abs() < 1e-12, "sum={s}");
        // Symmetry around center
        for i in 0..k.len() / 2 {
            assert!((k[i] - k[k.len() - 1 - i]).abs() < 1e-15);
        }
    }

    #[test]
    fn adaptive_filter_size_shrinks_for_tight_sigma() {
        let full = make_gaussian_kernel(50, 1.0);
        let fs = adaptive_filter_size(&full, GATK_BAND_PASS_MIN_PROB_TO_KEEP_IN_FILTER);
        assert_eq!(fs, 4);
        let k = make_gaussian_kernel(fs, 1.0);
        assert_eq!(k.len(), 9);
        assert!((k.iter().sum::<f64>() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn activity_state_offset_matches_gatk_semantics() {
        let st = ActivityProfileState::new("chr1", 1005, 0.5);
        assert_eq!(st.offset_from_region_start(1000), 5);
    }

    #[test]
    fn positive_sigma_rejects_non_positive() {
        assert!(PositiveSigma::try_new(0.0).is_none());
        assert!(PositiveSigma::try_new(-1.0).is_none());
        assert!(PositiveSigma::try_new(f64::NAN).is_none());
        assert_eq!(PositiveSigma::try_new(17.0).unwrap().get(), 17.0);
    }

    #[test]
    fn hq_soft_clip_evidence_carries_length_plain_does_not() {
        let plain = ActivityProfileState::new("chr1", 1, 0.1);
        assert_eq!(plain.evidence, ActivityEvidence::Plain);
        assert_eq!(plain.hq_soft_clip_bases(), None);
        let hq = ActivityProfileState::high_quality_soft_clips("chr1", 1, 0.1, 12);
        assert_eq!(
            hq.evidence,
            ActivityEvidence::HighQualitySoftClips { clip_bases: 12 }
        );
        assert_eq!(hq.hq_soft_clip_bases(), Some(12));
    }

    #[test]
    fn is_active_respects_strict_greater_than_threshold() {
        let st = ActivityProfileState::new("chr1", 1, 0.002);
        assert!(!st.is_active(0.002));
        assert!(st.is_active(0.001999));
    }

    #[test]
    fn hq_expansion_span_matches_max_propagation_plus_filter_cap() {
        let just = ActivityProfileState::high_quality_soft_clips("chr1", 5000, 0.5, 10_000);
        let max_r = 54u32;
        let v = activity_profile_base_process_state(&just, max_r, 20_000);
        assert_eq!(v.len(), 2 * (max_r as usize) + 1);
        assert_eq!(v[0].pos, 5000u64 - max_r as u64);
    }

    #[test]
    fn single_centered_spike_drops_left_tail_before_region_start() {
        // `ActivityProfile#incorporateSingleState` ignores offsets `< 0` vs `regionStartLoc`, so the
        // left half of the Gaussian is discarded when the first input sits at the profile origin.
        let params = BandPassActivityProfileParams {
            max_prob_propagation_distance: 50,
            active_prob_threshold: 0.002,
            max_filter_size: 1,
            sigma: PositiveSigma::try_new(1.0).unwrap(),
            adaptive_filter_size: false,
        };
        assert_eq!(params.resolved_filter_size(), 1);
        let mut prof = BandPassActivityProfile::new("chr1", 10_000, params);
        prof.add(ActivityProfileState::new("chr1", 5000, 1.0))
            .unwrap();
        let t = prof.total_smoothed_mass();
        // kernel half-width 1 at σ=1: retained mass = w_center + w_right = 1 - w_left
        let expected = 0.725931380938803;
        assert!((t - expected).abs() < 1e-12, "t={t}");
        assert_eq!(prof.len(), 2);
        assert_eq!(prof.states()[0].pos, 5000);
        assert_eq!(prof.states()[1].pos, 5001);
    }

    #[test]
    fn two_adjacent_spikes_accumulate_overlapping_windows() {
        let params = BandPassActivityProfileParams {
            max_prob_propagation_distance: 50,
            active_prob_threshold: 0.002,
            max_filter_size: 1,
            sigma: PositiveSigma::try_new(1.0).unwrap(),
            adaptive_filter_size: false,
        };
        let mut prof = BandPassActivityProfile::new("chr1", 10_000, params);
        prof.add(ActivityProfileState::new("chr1", 5000, 1.0))
            .unwrap();
        prof.add(ActivityProfileState::new("chr1", 5001, 1.0))
            .unwrap();
        assert_eq!(prof.len(), 3);
        let expected_total = 1.725931380938803;
        let t = prof.total_smoothed_mass();
        assert!(
            (t - expected_total).abs() < 1e-12,
            "unexpected total mass {t}"
        );
    }

    #[test]
    fn add_rejects_non_contiguous_locus() {
        let params = BandPassActivityProfileParams {
            max_prob_propagation_distance: 50,
            active_prob_threshold: 0.002,
            max_filter_size: 1,
            sigma: PositiveSigma::try_new(1.0).unwrap(),
            adaptive_filter_size: false,
        };
        let mut prof = BandPassActivityProfile::new("chr1", 10_000, params);
        prof.add(ActivityProfileState::new("chr1", 100, 1.0))
            .unwrap();
        let err = prof
            .add(ActivityProfileState::new("chr1", 102, 1.0))
            .unwrap_err();
        let msg = format!("{err}");
        assert!(msg.contains("not immediately after"));
    }

    #[test]
    fn pop_ready_regions_waits_until_safe_without_force_conversion() {
        let params = BandPassActivityProfileParams {
            max_prob_propagation_distance: 2,
            active_prob_threshold: 0.2,
            max_filter_size: 0,
            sigma: PositiveSigma::try_new(1.0).unwrap(),
            adaptive_filter_size: false,
        };
        let mut p = BandPassActivityProfile::new("chr1", 10_000, params);
        for i in 1..=6 {
            p.add(ActivityProfileState::new("chr1", i, 0.9)).unwrap();
        }
        // Need at least maxRegion + maxProp (= 5 + 2) states to finalize without force.
        let popped = p.pop_ready_regions(100, 2, 5, false).unwrap();
        assert!(popped.is_empty());
    }

    #[test]
    fn pop_ready_regions_cuts_active_region_at_local_minimum_when_maxed() {
        let params = BandPassActivityProfileParams {
            max_prob_propagation_distance: 0,
            active_prob_threshold: 0.2,
            max_filter_size: 0,
            sigma: PositiveSigma::try_new(1.0).unwrap(),
            adaptive_filter_size: false,
        };
        let mut p = BandPassActivityProfile::new("chr1", 10_000, params);
        // First 5 loci remain active (> threshold), forcing maxRegion cut logic.
        let probs = [0.9, 0.8, 0.7, 0.25, 0.6, 0.1];
        for (idx, &pr) in probs.iter().enumerate() {
            p.add(ActivityProfileState::new("chr1", (idx + 1) as u64, pr))
                .unwrap();
        }
        let popped = p.pop_ready_regions(2, 2, 5, true).unwrap();
        assert!(!popped.is_empty());
        assert_eq!(popped[0].start, 1);
        // maxRegion=5 and local minimum at locus 4 should become the cut.
        assert_eq!(popped[0].end, 4);
        assert!(popped[0].is_active);
    }

    #[test]
    fn pop_ready_regions_splits_active_then_inactive_and_updates_origin() {
        let params = BandPassActivityProfileParams {
            max_prob_propagation_distance: 0,
            active_prob_threshold: 0.2,
            max_filter_size: 0,
            sigma: PositiveSigma::try_new(1.0).unwrap(),
            adaptive_filter_size: false,
        };
        let mut p = BandPassActivityProfile::new("chr1", 10_000, params);
        let probs = [0.9, 0.8, 0.0, 0.0];
        for (idx, &pr) in probs.iter().enumerate() {
            p.add(ActivityProfileState::new("chr1", (idx + 1) as u64, pr))
                .unwrap();
        }
        let popped = p.pop_ready_regions(1, 1, 10, true).unwrap();
        assert_eq!(popped.len(), 2);
        assert_eq!(
            (popped[0].start, popped[0].end, popped[0].is_active),
            (1, 2, true)
        );
        assert_eq!(
            (popped[1].start, popped[1].end, popped[1].is_active),
            (3, 4, false)
        );
        assert!(p.is_empty());
        assert_eq!(p.region_start(), None);
    }

    #[test]
    fn padded_bounds_are_clipped_to_contig_edges() {
        let params = BandPassActivityProfileParams {
            max_prob_propagation_distance: 0,
            active_prob_threshold: 0.2,
            max_filter_size: 0,
            sigma: PositiveSigma::try_new(1.0).unwrap(),
            adaptive_filter_size: false,
        };
        let mut p = BandPassActivityProfile::new("chr1", 100, params);
        p.add(ActivityProfileState::new("chr1", 1, 0.9)).unwrap();
        let popped = p.pop_ready_regions(25, 1, 10, true).unwrap();
        assert_eq!(popped.len(), 1);
        assert_eq!(popped[0].padded_start, 1);
        assert_eq!(popped[0].padded_end, 26);
    }
}
