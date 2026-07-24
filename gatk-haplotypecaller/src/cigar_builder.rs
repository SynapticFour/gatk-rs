//! GATK `CigarBuilder` parity (merge ops, trim end deletions).

use crate::cigar::{Cigar, CigarElement, CigarOperator};

/// Final CIGAR from [`CigarBuilder`] including end-deletion trim accounting.
/// # Invariants
/// `cigar` is the merged HTSJDK-ordered elements after [`CigarBuilder::make`].
/// Leading/trailing deletion trim counts match GATK end-trim behavior when enabled.
/// # Ownership
/// Owns built [`Cigar`] plus trim counters.
/// # Mutation
/// Immutable result of a single builder pass.
/// # Biological assumptions
/// Trailing/leading deletions at alignment ends are artifacts, not biological indels.
/// # Java equivalence
/// GATK `CigarBuilder` output (`org.broadinstitute.hellbender.utils.sam.CigarBuilder`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CigarBuilderResult {
    pub cigar: Cigar,
    pub leading_deletions_removed: usize,
    pub trailing_deletions_removed: usize,
}

#[derive(PartialEq, Eq)]
enum Section {
    LeftHardClip,
    LeftSoftClip,
    Middle,
    RightSoftClip,
    RightHardClip,
}

/// Incremental CIGAR constructor with optional end-deletion trimming.
/// # Invariants
/// Maintains HTSJDK clip section ordering (hard → soft → middle → soft → hard).
/// `remove_deletions_at_ends` mirrors GATK builder trim mode when true.
/// # Ownership
/// Owns in-progress element buffer; [`Self::make`] transfers elements into [`Cigar`].
/// # Mutation
/// [`Self::add`] mutates internal state; [`Self::make`] clears the builder for reuse.
/// # Biological assumptions
/// Same as [`Cigar`]: alignment notation, not variant normalization.
/// # Java equivalence
/// GATK `CigarBuilder`.
pub struct CigarBuilder {
    elements: Vec<CigarElement>,
    last_operator: Option<CigarOperator>,
    section: Section,
    remove_deletions_at_ends: bool,
    leading_deletions_removed: usize,
    trailing_deletions_removed: usize,
    trailing_deletions_removed_in_make: usize,
}

impl CigarBuilder {
    pub fn new(remove_deletions_at_ends: bool) -> Self {
        Self {
            elements: Vec::new(),
            last_operator: None,
            section: Section::LeftHardClip,
            remove_deletions_at_ends,
            leading_deletions_removed: 0,
            trailing_deletions_removed: 0,
            trailing_deletions_removed_in_make: 0,
        }
    }

    pub fn default_trim() -> Self {
        Self::new(true)
    }

    pub fn add(&mut self, element: CigarElement) -> &mut Self {
        if element.length == 0 {
            return self;
        }
        let operator = element.operator;
        if self.remove_deletions_at_ends && operator == CigarOperator::Deletion {
            if self.last_operator.is_none() || self.last_operator.is_some_and(|o| o.is_clipping()) {
                self.leading_deletions_removed += element.length;
                return self;
            }
            if self.last_operator == Some(CigarOperator::Insertion)
                && (self.elements.len() == 1
                    || self
                        .elements
                        .get(self.elements.len().saturating_sub(2))
                        .is_some_and(|e| e.operator.is_clipping()))
            {
                self.leading_deletions_removed += element.length;
                return self;
            }
        }
        self.advance_section(operator);
        if Some(operator) == self.last_operator {
            let n = self.elements.len() - 1;
            self.elements[n].length += element.length;
        } else if self.last_operator.is_none() {
            self.elements.push(element);
            self.last_operator = Some(operator);
        } else if operator.is_clipping() {
            if self.remove_deletions_at_ends
                && self
                    .last_operator
                    .is_some_and(|o| !o.consumes_read_bases() && !o.is_clipping())
            {
                self.trailing_deletions_removed +=
                    self.elements.last().map(|e| e.length).unwrap_or(0);
                if let Some(last) = self.elements.last_mut() {
                    *last = element;
                }
                self.last_operator = Some(operator);
            } else if self.remove_deletions_at_ends && self.last_two_deletion_insertion() {
                self.trailing_deletions_removed += self.elements[self.elements.len() - 2].length;
                let ins = self.elements.pop().unwrap();
                self.elements.pop();
                self.elements.push(ins);
                self.elements.push(element);
                self.last_operator = Some(operator);
            } else {
                self.elements.push(element);
                self.last_operator = Some(operator);
            }
        } else if operator == CigarOperator::Deletion
            && self.last_operator == Some(CigarOperator::Insertion)
        {
            let size = self.elements.len();
            if size > 1 && self.elements[size - 2].operator == CigarOperator::Deletion {
                self.elements[size - 2].length += element.length;
            } else {
                let ins = self.elements.pop().unwrap();
                self.elements.push(element);
                self.elements.push(ins);
            }
        } else {
            self.elements.push(element);
            self.last_operator = Some(operator);
        }
        self
    }

    fn last_two_deletion_insertion(&self) -> bool {
        self.last_operator == Some(CigarOperator::Insertion)
            && self.elements.len() > 1
            && self.elements[self.elements.len() - 2].operator == CigarOperator::Deletion
    }

    fn advance_section(&mut self, operator: CigarOperator) {
        match operator {
            CigarOperator::HardClip => {
                if matches!(
                    self.section,
                    Section::LeftSoftClip | Section::Middle | Section::RightSoftClip
                ) {
                    self.section = Section::RightHardClip;
                }
            }
            CigarOperator::SoftClip => {
                if matches!(self.section, Section::LeftHardClip) {
                    self.section = Section::LeftSoftClip;
                } else if self.section == Section::Middle {
                    self.section = Section::RightSoftClip;
                }
            }
            _ => {
                if matches!(self.section, Section::LeftHardClip | Section::LeftSoftClip) {
                    self.section = Section::Middle;
                }
            }
        }
    }

    pub fn make(&mut self) -> Cigar {
        self.trailing_deletions_removed_in_make = 0;
        if self.remove_deletions_at_ends && self.last_operator == Some(CigarOperator::Deletion) {
            self.trailing_deletions_removed_in_make =
                self.elements.last().map(|e| e.length).unwrap_or(0);
            self.elements.pop();
        } else if self.remove_deletions_at_ends && self.last_two_deletion_insertion() {
            self.trailing_deletions_removed_in_make = self.elements[self.elements.len() - 2].length;
            self.elements.remove(self.elements.len() - 2);
        }
        Cigar {
            elements: std::mem::take(&mut self.elements),
        }
    }

    pub fn make_and_record(&mut self) -> CigarBuilderResult {
        let cigar = self.make();
        CigarBuilderResult {
            cigar,
            leading_deletions_removed: self.leading_deletions_removed,
            trailing_deletions_removed: self.trailing_deletions_removed
                + self.trailing_deletions_removed_in_make,
        }
    }
}
