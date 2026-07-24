# Phase C.1 — raw `ActivityProfileState` per locus

**Prerequisite:** Phase B complete — `./scripts/parity/run_hc_full_parity_phase_b.sh` must pass before relying on C gates.

L1 gate: `scripts/parity/run_hc_full_parity_c1_raw_activity.sh`

Producer: `hc_full_parity_gate_dump raw-activity <ref.fa> <bam|sam> <interval_cli>`

TSV columns: `contig`, `pos` (1-based), `active_prob`, `original_active_prob`, `kind` (`none` | `hq_soft_clips`).

Values are **pre–band-pass** pileup scores (same path as `evaluate_hc_activity_state` / GATK `isActive` input). L2 Java per-locus dump is optional (see `../java_dumps/README.md`).
