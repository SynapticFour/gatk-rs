//! HTSJDK-style CIGAR elements for assembly haplotype alignment.

/// SAM CIGAR operator subset used for haplotype alignment (M/I/D/S/H).
/// # Invariants
/// Consumption rules follow SAM: Match/Del consume ref; Match/Ins/SoftClip consume read.
/// # Ownership
/// [`Copy`] enum.
/// # Mutation
/// Immutable discriminant on [`CigarElement`].
/// # Biological assumptions
/// Encodes gapped alignment ops between haplotype/read and reference.
/// # Java equivalence
/// HTSJDK / GATK `CigarOperator` (assembly subset).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CigarOperator {
    Match,
    Insertion,
    Deletion,
    SoftClip,
    HardClip,
}

impl CigarOperator {
    pub fn consumes_reference_bases(self) -> bool {
        matches!(self, Self::Match | Self::Deletion)
    }

    pub fn consumes_read_bases(self) -> bool {
        matches!(self, Self::Match | Self::Insertion | Self::SoftClip)
    }

    pub fn is_clipping(self) -> bool {
        matches!(self, Self::SoftClip | Self::HardClip)
    }

    pub fn is_indel(self) -> bool {
        matches!(self, Self::Insertion | Self::Deletion)
    }

    pub fn is_alignment(self) -> bool {
        self == Self::Match
    }

    pub fn as_char(self) -> char {
        match self {
            Self::Match => 'M',
            Self::Insertion => 'I',
            Self::Deletion => 'D',
            Self::SoftClip => 'S',
            Self::HardClip => 'H',
        }
    }
}

/// One CIGAR operation with run length (HTSJDK-style subset used in assembly).
/// # Invariants
/// `length` > 0 in normalized [`Cigar`] sequences (builder may coalesce adjacent identical ops).
/// # Ownership
/// [`Copy`] pair of length + operator.
/// # Mutation
/// Immutable; [`Cigar::push`] may extend the last element instead of appending.
/// # Biological assumptions
/// Operators follow SAM spec consumption rules ([`CigarOperator::consumes_reference_bases`], etc.).
/// # Java equivalence
/// HTSJDK / GATK `CigarElement` (`org.broadinstitute.hellbender.utils.sam.CigarElement`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CigarElement {
    pub length: usize,
    pub operator: CigarOperator,
}

/// Ordered CIGAR run-length encoding for haplotype vs reference alignment.
/// # Invariants
/// Adjacent elements with the same operator are merged by [`Self::push`].
/// [`Self::reference_length`] / [`Self::read_length`] sum operator consumption over elements.
/// # Ownership
/// Owns `elements` vector; cheap to clone for haplotype snapshots.
/// # Mutation
/// Append-only via [`Self::push`]; callers may replace `elements` wholesale.
/// # Biological assumptions
/// Represents gapped alignment between haplotype read coordinates and reference span.
/// # Java equivalence
/// HTSJDK / GATK `Cigar` (`CigarUtils`, assembly sanity filters).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Cigar {
    pub elements: Vec<CigarElement>,
}

impl Cigar {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn push(&mut self, length: usize, operator: CigarOperator) {
        if length == 0 {
            return;
        }
        if let Some(last) = self.elements.last_mut() {
            if last.operator == operator {
                last.length += length;
                return;
            }
        }
        self.elements.push(CigarElement { length, operator });
    }

    pub fn reference_length(&self) -> usize {
        self.elements
            .iter()
            .map(|e| {
                if e.operator.consumes_reference_bases() {
                    e.length
                } else {
                    0
                }
            })
            .sum()
    }

    /// Longest insertion or deletion element (GATK assembly sanity filter).
    pub fn max_indel_element_length(&self) -> usize {
        self.elements
            .iter()
            .filter(|e| e.operator.is_indel())
            .map(|e| e.length)
            .max()
            .unwrap_or(0)
    }

    pub fn read_length(&self) -> usize {
        self.elements
            .iter()
            .map(|e| {
                if e.operator.consumes_read_bases() {
                    e.length
                } else {
                    0
                }
            })
            .sum()
    }

    pub fn to_gatk_string(&self) -> String {
        let mut out = String::new();
        for e in &self.elements {
            out.push_str(&format!("{}{}", e.length, e.operator.as_char()));
        }
        out
    }
}

pub fn length_on_reference(op: CigarOperator, len: usize) -> usize {
    if op.consumes_reference_bases() {
        len
    } else {
        0
    }
}

pub fn length_on_read(op: CigarOperator, len: usize) -> usize {
    if op.consumes_read_bases() {
        len
    } else {
        0
    }
}
