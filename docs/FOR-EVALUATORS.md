# For evaluators

What this tree claims, and what it does not. The only product-claim authority is [CLAIM_MATRIX.md](CLAIM_MATRIX.md). There is no gatk-rs claim in the website registry, so this page does not add a public status.

## Workspace version

`Cargo.toml` sets the workspace version to 0.1.0. That is the package version. It is not a public status label.

## What the default command runs

The default `HaplotypeCaller` command runs the `assembly-region-v1` path (`CallRegionArgs::strict_java()`): activity regions, local reassembly, PairHMM, genotyping, and VCF or gVCF emission. `CombineGVCFs`, `GenotypeGVCFs`, and `VariantFiltration` are real subcommands.

`gatk-haplotypecaller/src/run.rs` names that default a production default, and the VCF header records `##GATK_RS_HC_PIPELINE=assembly-region-v1`. That word is the name of the code path. It is not a clinical deployment and not a statement that the caller is validated for care.

## Yes rows

Product claims are only the **Yes** rows in [CLAIM_MATRIX.md](CLAIM_MATRIX.md). In short, those rows cover:

- the canonical mid-B HaplotypeCaller path on one ActiveFull region, `2:92317000-92319000` (three oracle SNPs);
- P12 L3/L4/L5 on `chr2:92300000-92350000` (66 Java sites);
- L2 synthetic gates, recorded **223/223 strict** on **2026-07-24** (that date is the matrix evidence; this page does not re-certify a later commit);
- a synthetic joint-genotype ladder, recommended **≤ 100** samples on a 10 kb interval;
- hard-filter FILTER parity on a synthetic SNP boundary set.

HG001 GIAB **ci-subset** F1 is a historical record at `f70abf86` (#155), pending revalidation. It is not a current-main **Yes** row. Full-autosome GIAB is not complete.

A CI badge is a workflow link, not those claims.

Waivers that limit an otherwise **Yes** row (W-H1, W-H3, W-L7-FORMAT, and the W-J* tables) are named in the claim matrix. Do not drop them when citing a gate.

## What is not GATK equivalence

| In this tree | Not in this tree |
|--------------|------------------|
| Germline HaplotypeCaller, then CombineGVCFs, GenotypeGVCFs, and hard `VariantFiltration` | BQSR, VQSR, Mutect2, gCNV/SV, Funcotator |
| The scoped gates above | Genome-wide or full-autosome HaplotypeCaller equivalence |
| One-sample HaplotypeCaller, then Combine | Java-style joint HaplotypeCaller that merges multiple `-I` BAMs into one read set |
| Hard filters (`--preset snp` / `indel`) | VQSR |
| Omitted `--pair-hmm` selects `FASTEST_AVAILABLE` (host SIMD when the build can use it, otherwise Logless scalar). Log10 scalar is an explicit `--pair-hmm LOG10_PAIRHMM` | A signed claim that the omitted flag matches GATK's log10 PairHMM, or a signed SIMD speed claim |

Also not claimed: bitwise-identical QUAL/FORMAT on every site, bamout, DRAGSTR, DRAGEN mode, allele-specific `AS_*` annotations, or `--assembly-region-out`. Multi-sample Combine loads gVCFs in memory and has no GenomicsDB path.

On the signed scopes, the caller is exercised against GATK 4.4.0.0 (source SHA `2dbc025821bc5f686c423ff332a41e6cef892a77`, image `us.gcr.io/broad-gatk/gatk:4.4.0.0`). Outside those scopes, agreement is unproven. The CLI accepts familiar `-R` / `-I` / `-O` / `-L` flags. `--java-options` is accepted and does not start a JVM.

Whole-codebase HaplotypeCaller parity is not established. The 6R/7R notes under [parity/](parity/) are investigation records, not extra Yes rows.

This is an independent research reimplementation. It is not affiliated with, endorsed by, or supported by the Broad Institute. "GATK" is a trademark of the Broad Institute. See [NOTICE.md](../NOTICE.md). It is not a drop-in replacement for Broad Institute GATK and it is not clinically validated.

## Contact

[contact@synapticfour.com](mailto:contact@synapticfour.com)
