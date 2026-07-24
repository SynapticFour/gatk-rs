//! Typed biallelic PL→GL shapes from Java HC VCF (RN-1).
//! These are **evidence-class** templates (AD / cluster allele geometry), not locus pins.
//! Absolute P12 sites remain characterization fixtures elsewhere (W-H1).

/// Observed Java HC biallelic PL vectors and their GL inverses (evidence-class templates).
/// # Invariants
/// Each variant's [`Self::pl`] / [`Self::gl`] tables are fixed Java VCF inverses (RN-1), except `HomRefTrap`.
/// Not locus pins — classify AD/PL geometry for sparse rescue paths.
/// # Ownership
/// [`Copy`] enum; GL/PL accessors return stack arrays or owned `Vec` via helpers.
/// # Mutation
/// Immutable template selection via [`Self::from_ad_pair`] / [`Self::from_pileup_depths`].
/// # Biological assumptions
/// Biallelic diploid site with informative read AD separating hom-alt vs het rescue shapes.
/// # Java equivalence
/// Rust-native evidence templates measured from Java HC VCF (L8 generalization); not a Java enum.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SparsePlShape {
    /// PL `90,6,0` → GL `[-9.0, -0.6, 0.0]` (hom-alt; QUAL≈78.32).
    HomAltStrong,
    /// PL `45,3,0` → GL `[-4.5, -0.3, 0.0]` (weak hom-alt / coupled cluster).
    HomAltWeak,
    /// PL `81,0,36` → GL `[-8.1, 0.0, -3.6]` (het).
    Het,
    /// PL `39,0,39` → GL `[-3.9, 0.0, -3.9]` (balanced het / CTC).
    HetBalanced,
    /// Synthetic fallback when no alt support (not a Java VCF rescue shape).
    HomRefTrap,
}

impl SparsePlShape {
    pub const fn gl(self) -> [f64; 3] {
        match self {
            Self::HomAltStrong => [-9.0, -0.6, 0.0],
            Self::HomAltWeak => [-4.5, -0.3, 0.0],
            Self::Het => [-8.1, 0.0, -3.6],
            Self::HetBalanced => [-3.9, 0.0, -3.9],
            Self::HomRefTrap => [-0.5, -2.0, -3.0],
        }
    }

    pub const fn pl(self) -> [i32; 3] {
        match self {
            Self::HomAltStrong => [90, 6, 0],
            Self::HomAltWeak => [45, 3, 0],
            Self::Het => [81, 0, 36],
            Self::HetBalanced => [39, 0, 39],
            Self::HomRefTrap => [0, 15, 25],
        }
    }

    pub fn gl_vec(self) -> Vec<f64> {
        self.gl().to_vec()
    }

    /// Java sparse rescue from informative AD pair (hom-alt vs het), not pileup depth alone.
    pub fn from_ad_pair(read_ref_ad: i32, read_alt_ad: i32) -> Option<Self> {
        let ra = read_alt_ad.max(0);
        let rr = read_ref_ad.max(0);
        if ra < 1 {
            return None;
        }
        if rr == 0 && ra >= 2 {
            Some(Self::HomAltStrong)
        } else if rr == 0 && ra == 1 {
            Some(Self::HomAltWeak)
        } else {
            Some(Self::Het)
        }
    }

    /// Pileup genotype when PairHMM is absent: alt-dominated → strong hom-alt, else het.
    pub fn from_pileup_depths(read_ref_ad: i32, read_alt_ad: i32) -> Self {
        let ref_ad = read_ref_ad.max(0);
        let alt_ad = read_alt_ad.max(0);
        let alt_dominated = alt_ad >= 2 && (ref_ad == 0 || alt_ad >= ref_ad.saturating_mul(4));
        if alt_dominated {
            Self::HomAltStrong
        } else if alt_ad >= 1 {
            Self::Het
        } else {
            Self::HomRefTrap
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ad_pair_maps_java_rescue_shapes() {
        assert_eq!(
            SparsePlShape::from_ad_pair(0, 2),
            Some(SparsePlShape::HomAltStrong)
        );
        assert_eq!(
            SparsePlShape::from_ad_pair(0, 1),
            Some(SparsePlShape::HomAltWeak)
        );
        assert_eq!(SparsePlShape::from_ad_pair(1, 2), Some(SparsePlShape::Het));
        assert_eq!(SparsePlShape::from_ad_pair(3, 0), None);
    }

    #[test]
    fn pl_gl_tables_match_java_vcf_inverses() {
        assert_eq!(SparsePlShape::HomAltStrong.pl(), [90, 6, 0]);
        assert_eq!(SparsePlShape::HomAltStrong.gl(), [-9.0, -0.6, 0.0]);
        assert_eq!(SparsePlShape::Het.pl(), [81, 0, 36]);
        assert_eq!(SparsePlShape::Het.gl(), [-8.1, 0.0, -3.6]);
        assert_eq!(SparsePlShape::HetBalanced.pl(), [39, 0, 39]);
    }
}
