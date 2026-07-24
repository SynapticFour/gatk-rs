# Phase D.2 — downsampling

Gate: `scripts/parity/run_hc_full_parity_d2_downsample.sh`

- **Positional:** `downsample-positional <sam> <cap> [non-random|random]` — pro `SAMRecord.getAlignmentStart()` (wie Java-Gate); Reads ohne zugewiesene Position (`ReadUtils.readHasNoAssignedPosition`) werden **nicht** reservoir-gekappt; optionaler Modus spiegelt `PositionalDownsampler(..., nonRandom)` wider.
- **Allele-biased (scaffold):** `downsample-allele <ref.fa> <sam> <contig> <pos1> <cap>` — round-robin ref/alt classes at one locus.

`parity/fixtures/hc-full-parity/d2/cases.tsv` Spalte `positional_mode` (9.) nur für `positional_*`-Zeilen; leer = non-random.
