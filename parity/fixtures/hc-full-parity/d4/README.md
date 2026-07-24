# Phase D.4 — read shard pipeline (IUPAC + filter + post-transform)

Producer: `hc_full_parity_gate_dump read-shard-pipeline <sam|bam>`

Mirrors GATK `MultiIntervalLocalReadShard` order **without** downsampler (D.2):

1. `HaplotypeCallerEngine.makeStandardHCReadTransformer()` — strict `IUPACReadTransformer`
2. `HaplotypeCallerEngine.makeStandardHCReadFilters()`
3. `ReadTransformer.identity()` post-filter hook (HC default)

Schema: `qname`, `flags`, `mapq`, `seq_raw`, `seq_after_pre`, `passes_hc_filter`, `seq_after_post`.

Include `@RG` with `SM:` on synthetic SAM inputs (WellformedReadFilter).
