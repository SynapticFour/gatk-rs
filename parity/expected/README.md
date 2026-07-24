# Golden / expected outputs (Phase 0 / Step 3)

This directory holds **committed** expected artifacts for differential tests (VCF/BAM/SAM snippets) once a tool reaches **bitwise** or **normalized** parity gates.

## Policy

- Filenames: `{tool}-{fixture}-{mode}.expected.vcf` (or `.sam`, `.bam` for small slices).
- Prefer **small** slices (tens of kb / few variants) to keep the repo lean.
- Document the exact Java command line used to produce each golden in the same PR.

## Current state

Smoke parity compares **live** Java vs Rust runs (`parity/reports/`). This directory also holds **small locked snippets** used as documentation / optional manual diff anchors:

- `countbases-chr1-1-16.lines.txt` — expected `CountBasesInReference -L chr1:1-16` histogram lines on `parity/fixtures/reference.fa` (Phase 1 / Step 24 anchor).
- `countbases-interval-list-full.lines.txt` — expected histogram lines for `-L parity/fixtures/regions.interval_list`.
- `countbases-whole-reference.lines.txt` — expected histogram lines for whole-reference traversal on `parity/fixtures/reference.fa`.
- `sample.strict.vcf` — strict VCF golden anchor used by `compare_vcf_strict.py` self-check.

VCF/BAM bitwise goldens for HC remain **future work** once genotyping output stabilizes.
