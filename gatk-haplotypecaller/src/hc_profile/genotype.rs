//! Genotyping / AD workload counters for production profiling.

#[derive(Debug, Clone, Default)]
pub struct GenotypeSiteSample {
    pub candidate_alleles: u64,
    pub genotype_states: u64,
    pub pl_vector_len: u64,
    /// Diploid sample count represented in this site call (HC single-sample → 1).
    pub samples: u64,
    pub wall_ns: u64,
    pub ad_wall_ns: u64,
    pub event_rebuild_wall_ns: u64,
    pub allele_map_wall_ns: u64,
    pub marginalize_wall_ns: u64,
    /// Diploid PL / genotype-index math (typically 3 states — not the dense hotspot).
    pub genotype_enum_wall_ns: u64,
}

#[derive(Debug, Default)]
pub struct GenotypeAgg {
    pub sites: u64,
    pub candidate_alleles_sum: u64,
    pub genotype_states_sum: u64,
    pub pl_vector_len_sum: u64,
    pub samples_sum: u64,
    pub wall_ns: u64,
    pub ad_wall_ns: u64,
    pub event_rebuild_wall_ns: u64,
    pub allele_map_wall_ns: u64,
    pub marginalize_wall_ns: u64,
    pub genotype_enum_wall_ns: u64,
}

impl GenotypeAgg {
    pub fn add(&mut self, s: GenotypeSiteSample) {
        self.sites += 1;
        self.candidate_alleles_sum += s.candidate_alleles;
        self.genotype_states_sum += s.genotype_states;
        self.pl_vector_len_sum += s.pl_vector_len;
        self.samples_sum += s.samples;
        self.wall_ns += s.wall_ns;
        self.ad_wall_ns += s.ad_wall_ns;
        self.event_rebuild_wall_ns += s.event_rebuild_wall_ns;
        self.allele_map_wall_ns += s.allele_map_wall_ns;
        self.marginalize_wall_ns += s.marginalize_wall_ns;
        self.genotype_enum_wall_ns += s.genotype_enum_wall_ns;
    }

    pub fn time_per_site_ns(&self) -> f64 {
        if self.sites == 0 {
            0.0
        } else {
            self.wall_ns as f64 / self.sites as f64
        }
    }

    pub fn time_per_state_ns(&self) -> f64 {
        if self.genotype_states_sum == 0 {
            0.0
        } else {
            self.wall_ns as f64 / self.genotype_states_sum as f64
        }
    }
}
