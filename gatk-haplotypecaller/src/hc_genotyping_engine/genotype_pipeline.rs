/// L14-B: thin production genotyping pipeline surface (map → … → filter).
/// Stages are owned by dedicated types; this trait documents the composition for
/// tests and future alternate engines. Default production wiring uses
/// [`StrictJavaGenotypingPipeline`].
pub(crate) trait GenotypingPipeline {
    /// Named stage order for production `enable_java_strict` sites.
    fn stage_names(&self) -> &'static [&'static str] {
        &[
            "map",
            "early_template",
            "pileup_rescue",
            "score",
            "reshape",
            "finalize",
            "filter",
        ]
    }
}

/// Production pipeline: [`SiteMap`] → [`SiteEarlyTemplate`] → [`SitePileupRescue`]
/// → [`SiteScore`] → [`SiteReshape`] → [`GenotypeFinalize`] → emit filter.
pub(crate) struct StrictJavaGenotypingPipeline;

impl GenotypingPipeline for StrictJavaGenotypingPipeline {}
