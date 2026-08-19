//! Named production stages for HC profiling.

use std::time::Duration;

/// Coarse production stages (stable JSON keys).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Stage {
    InputBamDecode,
    ReadPreprocess,
    ActiveRegionConstruction,
    EventDiscovery,
    AssemblyGraphConstruction,
    GraphPruning,
    HaplotypeGeneration,
    SmithWaterman,
    PairHmm,
    LikelihoodProcessing,
    GenotypeAssignment,
    AdAnnotation,
    VcfEmission,
    Synchronization,
    Allocations,
    /// Catch-all for TRACE phases that do not map cleanly.
    Other,
}

impl Stage {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::InputBamDecode => "input_bam_decode",
            Self::ReadPreprocess => "read_preprocessing",
            Self::ActiveRegionConstruction => "active_region_construction",
            Self::EventDiscovery => "event_discovery",
            Self::AssemblyGraphConstruction => "assembly_graph_construction",
            Self::GraphPruning => "graph_pruning",
            Self::HaplotypeGeneration => "haplotype_generation",
            Self::SmithWaterman => "smith_waterman",
            Self::PairHmm => "pairhmm",
            Self::LikelihoodProcessing => "likelihood_processing",
            Self::GenotypeAssignment => "genotype_assignment",
            Self::AdAnnotation => "ad_annotation",
            Self::VcfEmission => "vcf_emission",
            Self::Synchronization => "synchronization_waiting",
            Self::Allocations => "allocations",
            Self::Other => "other",
        }
    }

    pub fn all() -> &'static [Stage] {
        &[
            Self::InputBamDecode,
            Self::ReadPreprocess,
            Self::ActiveRegionConstruction,
            Self::EventDiscovery,
            Self::AssemblyGraphConstruction,
            Self::GraphPruning,
            Self::HaplotypeGeneration,
            Self::SmithWaterman,
            Self::PairHmm,
            Self::LikelihoodProcessing,
            Self::GenotypeAssignment,
            Self::AdAnnotation,
            Self::VcfEmission,
            Self::Synchronization,
            Self::Allocations,
            Self::Other,
        ]
    }

    /// Map `HC_RSS_TRACE phase=` names onto coarse stages.
    pub fn from_trace_phase(phase: &str) -> Option<Stage> {
        let p = phase;
        if p.starts_with("rt_build")
            || p.starts_with("rt_graph")
            || p.starts_with("rt_first")
            || p.starts_with("merge_rt")
            || p.contains("graph_build")
        {
            return Some(Self::AssemblyGraphConstruction);
        }
        if p.starts_with("seq_") || p.contains("prune") || p.contains("simplify") {
            return Some(Self::GraphPruning);
        }
        if p.starts_with("kbest") || p.contains("extract") {
            return Some(Self::HaplotypeGeneration);
        }
        if p.starts_with("rt_dangling") || p.contains("dangling") {
            return Some(Self::SmithWaterman);
        }
        if p.contains("pairhmm") || p == "before_pairhmm" || p == "after_pairhmm" {
            return Some(Self::PairHmm);
        }
        if p.contains("realign") || p.contains("sw_") {
            return Some(Self::SmithWaterman);
        }
        if p.contains("genotype") || p.contains("assign") {
            return Some(Self::GenotypeAssignment);
        }
        if p.contains("spine")
            || p.contains("event")
            || p.contains("emap")
            || p.contains("variation")
            || p.starts_with("prep_trim")
            || p.starts_with("prep_post")
            || p.starts_with("prep_parity")
            || p.starts_with("prep_early_allele")
        {
            return Some(Self::EventDiscovery);
        }
        if p.contains("emit") || p.contains("vcf") || p.contains("gvcf") {
            return Some(Self::VcfEmission);
        }
        if p.contains("finalize") || p == "before_finalize" || p == "after_finalize" {
            return Some(Self::ReadPreprocess);
        }
        if p == "after_assemble" || p.starts_with("assemble") {
            return Some(Self::AssemblyGraphConstruction);
        }
        if p == "run_start" {
            return None;
        }
        Some(Self::Other)
    }
}

#[derive(Debug, Default, Clone)]
pub struct StageStats {
    pub calls: u64,
    pub wall_ns: u64,
    pub cpu_ns: u64,
    pub cpu_samples: u64,
    pub alloc_bytes: u64,
    pub alloc_events: u64,
}

impl StageStats {
    pub fn add(&mut self, wall: Duration, cpu: Option<Duration>) {
        self.calls += 1;
        self.wall_ns = self.wall_ns.saturating_add(wall.as_nanos() as u64);
        if let Some(c) = cpu {
            self.cpu_ns = self.cpu_ns.saturating_add(c.as_nanos() as u64);
            self.cpu_samples += 1;
        }
    }

    pub fn add_wall_only(&mut self, wall: Duration) {
        self.calls += 1;
        self.wall_ns = self.wall_ns.saturating_add(wall.as_nanos() as u64);
    }

    pub fn avg_wall_ns(&self) -> f64 {
        if self.calls == 0 {
            0.0
        } else {
            self.wall_ns as f64 / self.calls as f64
        }
    }
}

#[derive(Debug, Default)]
pub struct StageAgg {
    pub stats: [StageStats; 16],
}

impl StageAgg {
    fn idx(stage: Stage) -> usize {
        match stage {
            Stage::InputBamDecode => 0,
            Stage::ReadPreprocess => 1,
            Stage::ActiveRegionConstruction => 2,
            Stage::EventDiscovery => 3,
            Stage::AssemblyGraphConstruction => 4,
            Stage::GraphPruning => 5,
            Stage::HaplotypeGeneration => 6,
            Stage::SmithWaterman => 7,
            Stage::PairHmm => 8,
            Stage::LikelihoodProcessing => 9,
            Stage::GenotypeAssignment => 10,
            Stage::AdAnnotation => 11,
            Stage::VcfEmission => 12,
            Stage::Synchronization => 13,
            Stage::Allocations => 14,
            Stage::Other => 15,
        }
    }

    pub fn add(&mut self, stage: Stage, wall: Duration, cpu: Option<Duration>) {
        self.stats[Self::idx(stage)].add(wall, cpu);
    }

    pub fn add_wall_only(&mut self, stage: Stage, wall: Duration) {
        self.stats[Self::idx(stage)].add_wall_only(wall);
    }

    pub fn add_alloc(&mut self, stage: Stage, bytes: u64, events: u64) {
        let s = &mut self.stats[Self::idx(stage)];
        s.alloc_bytes = s.alloc_bytes.saturating_add(bytes);
        s.alloc_events = s.alloc_events.saturating_add(events);
        // Also roll into Allocations bucket.
        let a = &mut self.stats[Self::idx(Stage::Allocations)];
        a.alloc_bytes = a.alloc_bytes.saturating_add(bytes);
        a.alloc_events = a.alloc_events.saturating_add(events);
    }

    pub fn get(&self, stage: Stage) -> &StageStats {
        &self.stats[Self::idx(stage)]
    }

    pub fn stats_mut(&mut self, stage: Stage) -> &mut StageStats {
        &mut self.stats[Self::idx(stage)]
    }
}
