use crate::parallel::ParallelConfig;

/// Async FASTA loader using Tokio filesystem I/O.
/// # Invariants
/// Parses minimal `>` header / sequence lines into name→sequence map.
/// # Ownership
/// Owns config; returns owned sequence vectors to caller.
/// # Mutation
/// Async methods borrow `&self`; no internal mutation.
/// # Biological assumptions
/// Standard FASTA text; concatenates wrapped sequence lines per record.
/// # Java equivalence
/// Similar role to htsjdk `ReferenceSequenceFile` async wrapper (Rust-native).
pub struct AsyncGenomicFileProcessor {
    _config: ParallelConfig,
}

impl AsyncGenomicFileProcessor {
    pub fn new(config: ParallelConfig) -> gatk_common::GatkResult<Self> {
        Ok(Self { _config: config })
    }

    pub async fn process_fasta_async(
        &self,
        path: &str,
    ) -> gatk_common::GatkResult<Vec<(String, Vec<u8>)>> {
        let content = tokio::fs::read_to_string(path)
            .await
            .map_err(|e| gatk_common::GatkError::io(format!("Failed to read FASTA: {path}"), e))?;

        let mut out: Vec<(String, Vec<u8>)> = Vec::new();
        let mut current_name: Option<String> = None;
        let mut current_seq: Vec<u8> = Vec::new();

        for line in content.lines() {
            if let Some(rest) = line.strip_prefix('>') {
                if let Some(name) = current_name.take() {
                    out.push((name, std::mem::take(&mut current_seq)));
                }
                current_name = Some(rest.trim().to_string());
            } else if !line.trim().is_empty() {
                current_seq.extend_from_slice(line.trim().as_bytes());
            }
        }

        if let Some(name) = current_name {
            out.push((name, current_seq));
        }

        Ok(out)
    }
}
