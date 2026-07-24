//! Parity v1 site annotations.

pub mod engine;
pub mod plugins;

pub use engine::{
    annotate_parity_v1_site, AnnotatedSite, VariantAnnotationContext, PARITY_V1_FORMAT_KEYS,
    PARITY_V1_INFO_KEYS,
};
