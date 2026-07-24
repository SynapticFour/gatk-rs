# Phase C.3 — per-locus binary isActive

**Prerequisite:** Phase B + C.1/C.2 recommended.

L1 gate: `scripts/parity/run_hc_full_parity_c3_active_locus.sh`

Producer: `hc_full_parity_gate_dump active-locus <ref.fa> <bam|sam> <interval_cli>`

TSV: `contig`, `pos`, `is_active` (`true` if smoothed `active_prob >` GATK `active-probability-threshold`, strict `>`).
