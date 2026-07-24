//! PairHMM likelihood dumps for L2 parity.

use crate::pairhmm::{pairhmm_log10_likelihood, PairHmmInput, PairHmmParams};
use gatk_common::{GatkError, GatkResult};
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

/// One PairHMM parity fixture case loaded from TSV (dumps).
/// # Invariants
/// `read_base_quals.len` matches `read_bases.len` when loaded from valid fixtures.
/// # Ownership
/// Owns case id, sequences, quals, and haplotype string.
/// # Mutation
/// Immutable after TSV load.
/// # Biological assumptions
/// None — deterministic test vector for likelihood parity.
/// # Java equivalence
/// Rust-native dump row; cases mirror GATK PairHMM unit/parity inputs.
#[derive(Debug, Clone)]
pub struct PairHmmCaseRow {
    pub case_id: String,
    pub read_bases: String,
    pub read_base_quals: Vec<u8>,
    pub read_mapq: u8,
    pub haplotype: String,
}

pub fn load_pairhmm_cases_tsv(path: &Path) -> GatkResult<Vec<PairHmmCaseRow>> {
    let f = File::open(path)
        .map_err(|e| GatkError::generic(format!("open {}: {e}", path.display())))?;
    let reader = BufReader::new(f);
    let mut out = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(|e| GatkError::generic(format!("read {}: {e}", path.display())))?;
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let mut parts = t.split('\t');
        let case_id = parts
            .next()
            .ok_or_else(|| GatkError::argument("pairhmm cases: missing case_id"))?
            .to_string();
        let read_bases = parts
            .next()
            .ok_or_else(|| GatkError::argument("pairhmm cases: missing read_bases"))?
            .to_string();
        let quals_str = parts
            .next()
            .ok_or_else(|| GatkError::argument("pairhmm cases: missing read_base_quals"))?;
        let read_base_quals: Vec<u8> = quals_str
            .split(',')
            .map(|q| {
                q.parse::<u8>()
                    .map_err(|_| GatkError::argument(format!("invalid qual {q}")))
            })
            .collect::<Result<Vec<_>, _>>()?;
        let read_mapq: u8 = parts
            .next()
            .ok_or_else(|| GatkError::argument("pairhmm cases: missing read_mapq"))?
            .parse()
            .map_err(|_| GatkError::argument("pairhmm cases: invalid read_mapq"))?;
        let haplotype = parts
            .next()
            .ok_or_else(|| GatkError::argument("pairhmm cases: missing haplotype"))?
            .to_string();
        out.push(PairHmmCaseRow {
            case_id,
            read_bases,
            read_base_quals,
            read_mapq,
            haplotype,
        });
    }
    Ok(out)
}

/// TSV: `case_id\tlog10_likelihood` (one row per case in fixture).
pub fn dump_pairhmm_likelihoods_tsv(cases_path: &Path, out: &mut impl Write) -> GatkResult<()> {
    let cases = load_pairhmm_cases_tsv(cases_path)?;
    let params = PairHmmParams::default();
    writeln!(out, "case_id\tlog10_likelihood")
        .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    for row in &cases {
        let ll = pairhmm_log10_likelihood(
            &PairHmmInput {
                read_bases: row.read_bases.clone(),
                read_base_quals: row.read_base_quals.clone(),
                read_mapping_quality: row.read_mapq,
                haplotype_bases: row.haplotype.clone(),
            },
            &params,
        )?;
        writeln!(out, "{}\t{}", row.case_id, ll)
            .map_err(|e| GatkError::generic(format!("write tsv: {e}")))?;
    }
    Ok(())
}
