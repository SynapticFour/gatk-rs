# Phase B.2.2 — empty locus / zero-depth pileup (L2)

Gate: `locus-pileup` rows with `pileup_depth=0` where no read covers the base, including a **read-free BAM** (header only).

Rust: `dump_locus_pileup_tsv` · Java: `HcFullParityGateDump locus-pileup`.

Frozen L2: `parity/fixtures/hc-full-parity/java_dumps/b2-empty-locus/<case_id>_<PIN>.tsv`.
