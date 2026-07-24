use crate::parallel::ParallelConfig;
use futures::future::try_join_all;

/// Concurrent file read/write helper built on Tokio.
/// # Invariants
/// One async task per path in batch APIs; failures fail entire batch.
/// # Ownership
/// Owns config; returns owned `(path, bytes)` tuples.
/// # Mutation
/// Async methods borrow `&self`.
/// # Biological assumptions
/// None (generic I/O); paths typically FASTA/BAM/VCF inputs.
/// # Java equivalence
/// None / Rust-native.
pub struct TokioProcessor {
    _config: ParallelConfig,
}

impl TokioProcessor {
    pub fn new(config: ParallelConfig) -> gatk_common::GatkResult<Self> {
        Ok(Self { _config: config })
    }

    pub async fn read_files_concurrent(
        &self,
        paths: Vec<String>,
    ) -> gatk_common::GatkResult<Vec<(String, Vec<u8>)>> {
        let tasks = paths.into_iter().map(|path| async move {
            let bytes = tokio::fs::read(&path).await.map_err(|e| {
                gatk_common::GatkError::io(format!("Failed to read file {path}"), e)
            })?;
            Ok::<_, gatk_common::GatkError>((path, bytes))
        });
        try_join_all(tasks).await
    }

    pub async fn write_files_concurrent(
        &self,
        files: Vec<(String, Vec<u8>)>,
    ) -> gatk_common::GatkResult<()> {
        let tasks = files.into_iter().map(|(path, bytes)| async move {
            tokio::fs::write(&path, bytes)
                .await
                .map_err(|e| gatk_common::GatkError::io(format!("Failed to write file {path}"), e))
        });
        try_join_all(tasks).await?;
        Ok(())
    }
}
