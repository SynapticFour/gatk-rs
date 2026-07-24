# Phase C.2 — band-pass smoothed activity

**Prerequisite:** `./scripts/parity/run_hc_full_parity_phase_b.sh`

L1 gate: `scripts/parity/run_hc_full_parity_c2_smoothed_activity.sh`

Producer: `hc_full_parity_gate_dump smoothed-activity <ref.fa> <bam|sam> <interval_cli>`

TSV: `contig`, `pos`, `smoothed_active_prob` (after `BandPassActivityProfile`, HC defaults).
