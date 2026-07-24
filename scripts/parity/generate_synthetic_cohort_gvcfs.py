#!/usr/bin/env python3
"""Generate a synthetic multi-sample gVCF cohort for CombineGVCFs scale tests.

Builds N single-sample gVCFs on a contig by expanding the checked-in
`parity/combine_gvcfs/mini/` pattern (hom-ref blocks + heterozygous SNPs with
`<NON_REF>`). Variant positions are deterministic from the sample index so
Java↔Rust parity stays reproducible.

Each SNP locus uses **one shared ALT** for all carriers (presence varies). That
keeps the scale ladder focused on N-sample merge cost / GT parity, without
conflating it with Java's last→first multi-allelic ALT ordering (covered by the
2-sample mini gate + unit tests).

This is **not** a substitute for GIAB WGS cohorts — it exercises Combine/Genotype
scaling with sample count on a fixed interval.

Usage:
  python3 scripts/parity/generate_synthetic_cohort_gvcfs.py \\
    --n-samples 50 \\
    --out-dir parity/combine_gvcfs/cohort_scale/n50
"""
from __future__ import annotations

import argparse
import json
import pathlib
import textwrap


HEADER = """\
##fileformat=VCFv4.2
##contig=<ID=chr1,length={contig_len}>
##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}
"""

# Compact fixture (mini-compatible). Used when --snp-count is small / default override.
MINI_SNP_LOCI = [
    (10, "A", "G"),
    (25, "C", "T"),
    (40, "G", "A"),
    (55, "T", "C"),
    (70, "A", "C"),
    (85, "C", "G"),
    (100, "G", "T"),
    (115, "T", "A"),
]


def build_snp_loci(contig_len: int, snp_count: int, snp_stride: int) -> list[tuple[int, str, str]]:
    """Return (pos, ref, alt) list. Shared ALT per locus across all samples."""
    if snp_count <= len(MINI_SNP_LOCI) and contig_len <= 160:
        return [(p, r, a) for p, r, a in MINI_SNP_LOCI if p < contig_len]
    alts = ("G", "T", "C", "A")
    loci: list[tuple[int, str, str]] = []
    pos = 10
    i = 0
    while len(loci) < snp_count and pos < contig_len:
        ref = "ACGT"[(pos - 1) % 4]
        alt = alts[i % len(alts)]
        if alt == ref:
            alt = alts[(i + 1) % len(alts)]
        loci.append((pos, ref, alt))
        pos += snp_stride
        i += 1
    return loci


def make_contig_bases(contig_len: int, snp_loci: list[tuple[int, str, str]]) -> str:
    bases = ["ACGT"[i % 4] for i in range(contig_len)]
    if contig_len >= 1:
        bases[0] = "A"
    for pos, ref, _alt in snp_loci:
        if 1 <= pos <= contig_len:
            bases[pos - 1] = ref
    return "".join(bases)


def write_ref_fa(path: pathlib.Path, bases: str) -> None:
    path.write_text(f">chr1\n{bases}\n", encoding="utf-8")


def write_ref_dict(path: pathlib.Path, contig_len: int) -> None:
    path.write_text(
        textwrap.dedent(
            f"""\
            @HD\tVN:1.6\tSO:unsorted
            @SQ\tSN:chr1\tLN:{contig_len}\tM5:00000000000000000000000000000000\tUR:file:ref.fa
            """
        ),
        encoding="utf-8",
    )


def diploid_pl_het(gq: int = 99) -> str:
    # GT 0/1 vs REF/ALT/<NON_REF> → 6 PL entries (diploid, 3 alleles).
    return f"{gq},0,{gq},{gq},{gq},{gq}"


def diploid_pl_homref(gq: int) -> str:
    return f"0,{gq},{gq * 10}"


def emit_sample_gvcf(
    path: pathlib.Path,
    sample: str,
    sample_idx: int,
    contig_len: int,
    bases: str,
    snp_loci: list[tuple[int, str, str]],
) -> None:
    # HEADER already contains newlines; strip before joining body lines so we do not
    # insert a blank record line (Java Tribble / Rust both reject that).
    lines = HEADER.format(contig_len=contig_len, sample=sample).rstrip("\n").split("\n")
    first_snp = snp_loci[0][0] if snp_loci else contig_len + 1
    if first_snp > 1:
        end = first_snp - 1
        lines.append(
            f"chr1\t1\t.\t{bases[0]}\t<NON_REF>\t.\t.\tEND={end}\t"
            f"GT:GQ:DP:MIN_DP:PL\t0/0:{80 + (sample_idx % 15)}:{20 + (sample_idx % 10)}:"
            f"{15 + (sample_idx % 8)}:{diploid_pl_homref(60 + (sample_idx % 20))}"
        )
        cursor = first_snp
    else:
        cursor = 1

    for pos, ref, alt in snp_loci:
        if pos >= contig_len:
            break
        if cursor < pos:
            end = pos - 1
            gq = 40 + (sample_idx % 30)
            lines.append(
                f"chr1\t{cursor}\t.\t{bases[cursor - 1]}\t<NON_REF>\t.\t.\tEND={end}\t"
                f"GT:GQ:DP:MIN_DP:PL\t0/0:{gq}:{12 + (sample_idx % 8)}:"
                f"{8 + (sample_idx % 5)}:{diploid_pl_homref(gq)}"
            )
        # ~70% of samples carry a het at this locus; rest stay hom-ref through the base.
        carries = (sample_idx * 7 + pos) % 10 < 7
        if carries:
            ad_ref = 8 + (sample_idx % 6)
            ad_alt = 8 + ((sample_idx * 3) % 6)
            dp = ad_ref + ad_alt
            lines.append(
                f"chr1\t{pos}\t.\t{ref}\t{alt},<NON_REF>\t.\t.\t.\t"
                f"GT:AD:DP:GQ:PL\t0/1:{ad_ref},{ad_alt},0:{dp}:99:{diploid_pl_het()}"
            )
        else:
            lines.append(
                f"chr1\t{pos}\t.\t{ref}\t<NON_REF>\t.\t.\tEND={pos}\t"
                f"GT:GQ:DP:MIN_DP:PL\t0/0:50:10:8:{diploid_pl_homref(50)}"
            )
        cursor = pos + 1

    if cursor <= contig_len:
        lines.append(
            f"chr1\t{cursor}\t.\t{bases[cursor - 1]}\t<NON_REF>\t.\t.\tEND={contig_len}\t"
            f"GT:GQ:DP:MIN_DP:PL\t0/0:40:10:8:{diploid_pl_homref(40)}"
        )

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n-samples", type=int, required=True, help="Cohort size (e.g. 50)")
    ap.add_argument("--out-dir", type=pathlib.Path, required=True)
    ap.add_argument(
        "--contig-len",
        type=int,
        default=10000,
        help="Contig length (default 10000 for measurable Combine wall vs N)",
    )
    ap.add_argument(
        "--snp-count",
        type=int,
        default=400,
        help="Number of SNP loci (default 400)",
    )
    ap.add_argument(
        "--snp-stride",
        type=int,
        default=20,
        help="Bases between SNP loci when generating the dense ladder",
    )
    args = ap.parse_args()
    if args.n_samples < 2:
        raise SystemExit("--n-samples must be >= 2")

    out = args.out_dir
    out.mkdir(parents=True, exist_ok=True)
    samples_dir = out / "samples"
    samples_dir.mkdir(exist_ok=True)

    snp_loci = build_snp_loci(args.contig_len, args.snp_count, args.snp_stride)
    bases = make_contig_bases(args.contig_len, snp_loci)
    write_ref_fa(out / "ref.fa", bases)
    write_ref_dict(out / "ref.dict", args.contig_len)
    (out / "ref.fa.fai").write_text(
        f"chr1\t{args.contig_len}\t6\t{args.contig_len}\t{args.contig_len + 1}\n"
    )

    sample_names: list[str] = []
    for i in range(args.n_samples):
        # First three mimic GIAB trio naming for dashboard readability; rest synthetic.
        if i == 0:
            name = "HG002"
        elif i == 1:
            name = "HG003"
        elif i == 2:
            name = "HG004"
        else:
            name = f"SYN{i:03d}"
        sample_names.append(name)
        emit_sample_gvcf(
            samples_dir / f"{name}.g.vcf", name, i, args.contig_len, bases, snp_loci
        )

    manifest = {
        "n_samples": args.n_samples,
        "samples": sample_names,
        "contig": "chr1",
        "contig_len": args.contig_len,
        "snp_loci": [p for p, _, _ in snp_loci],
        "snp_count": len(snp_loci),
        "shared_alt_per_locus": True,
        "notes": (
            "Synthetic cohort for CombineGVCFs→GenotypeGVCFs sample-count scaling. "
            "Not a GIAB truth set. One shared ALT per locus; carrier set varies by sample."
        ),
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    (out / "samples.list").write_text(
        "\n".join(str(samples_dir / f"{s}.g.vcf") for s in sample_names) + "\n"
    )
    print(f"[cohort-gen] wrote {args.n_samples} gVCFs ({len(snp_loci)} SNPs) → {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
