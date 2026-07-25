//! Java-compatibility / P12-interval semantics (waivers **W-H1**, **W-H3**).
//! This module holds predicates and band tables derived from GATK 4.4 observable
//! behavior on the NA12878 P12 validation slice (`chr2:92300000–92350000`).
//! **Not a claim of genome-wide GATK equivalence.** Prefer phenotype predicates over
//! new coordinate literals. See `docs/ARCHITECTURE.md` and
//! `docs/CLAIM_MATRIX.md`.
//! Sprint J: `coupled_indel` is the preferred recognition path for cluster indels;
//! absolute band tables in `java_hc_site_semantics` are oracle/waiver surfaces.

pub mod coupled_indel;
pub mod java_hc_site_semantics;

pub use coupled_indel::{
    is_coupled_indel_for_genotyping, is_coupled_indel_member, is_ctc_del_for_genotyping,
    CoupledIndelCluster, COUPLED_INDEL_PARTNER_OFFSET,
};
