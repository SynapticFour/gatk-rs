# B.2.1 — Locus pileup depth (`LocusIteratorByState` + `IntervalAlignmentContextIterator`)

## Artifact (`*.tsv`)

`hc_full_parity_gate_dump locus-pileup <ref.fa> <bam> <interval_cli> [padding]`

- Header: `contig\tpos\tpileup_depth`
- One row per 1-based position in user `-L` intervals (including empty loci with depth `0`).
- Reads are loaded per contig shard with GATK padding (default **100**); pileup depth matches Java `AlignmentContext.size()`.

## L1 / L2

- **L1:** Rust vs `expected/*.tsv` (`run_hc_full_parity_b2_locus.sh`).
- **L2:** Rust vs `parity/fixtures/hc-full-parity/java_dumps/b2-locus/*_<pin>.tsv`.
