#!/usr/bin/env python3
"""Generate g2_subset_live_quad4.{fa,sam} for live 4→2 hap-score trim (sensitive assembly).

Four SNP branches at 15/35/55/75 on an 86 bp contig. The v75 cluster uses a 12 bp alt
window plus 8 N padding (20M) so samtools can index the BAM while assembly still recovers
four haplotypes under the sensitive profile.
"""

from __future__ import annotations

from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
CONTIG = "chrQuad4"
REF = (
    "TGCATGACTGATCGTACGATTCGAGCTAGTCGATCGATGCTAGCTAGGCTAACGTTAGCTAGTAACTGATCGATCGATACGTACGT"
)
SNP_LOCI = [(15, 5, "A"), (35, 5, "C"), (55, 5, "G"), (75, 5, "T")]


def window(ref: str, pos: int, length: int) -> str:
    return ref[pos - 1 : pos - 1 + length]


def sam_line(qname: str, pos: int, cigar: str, seq: str) -> str:
    if cigar.endswith("M") and cigar[:-1].isdigit():
        assert int(cigar[:-1]) == len(seq), (qname, pos, cigar, len(seq))
    return (
        f"{qname}\t0\t{CONTIG}\t{pos}\t60\t{cigar}\t*\t0\t0\t{seq}\t"
        f"{'F' * len(seq)}\tRG:Z:rg1\n"
    )


def main() -> None:
    assert len(REF) == 86
    specs: list[tuple[str, int, str, str]] = []
    for i in range(8):
        specs.append((f"r{i}", 10, "30M", window(REF, 10, 30)))
    for pos, alt_idx, alt in SNP_LOCI:
        if pos == 75:
            alt_bases = list(window(REF, pos, 12))
            alt_bases[alt_idx] = alt
            seq = "".join(alt_bases) + ("N" * 8)
            cigar = "20M"
        else:
            alt_bases = list(window(REF, pos, 20))
            alt_bases[alt_idx] = alt
            seq = "".join(alt_bases)
            cigar = "20M"
        for i in range(5):
            specs.append((f"v{pos}_{i}", pos, cigar, seq))

    fa = REPO / "parity/fixtures/g2_subset_live_quad4.fa"
    sam = REPO / "parity/fixtures/g2_subset_live_quad4.sam"
    fa.write_text(f">{CONTIG}\n{REF}\n")
    hdr = (
        f"@HD\tVN:1.6\tSO:coordinate\n"
        f"@SQ\tSN:{CONTIG}\tLN:{len(REF)}\n"
        f"@RG\tID:rg1\tSM:s1\tPL:ILLUMINA\n"
    )
    sam.write_text(hdr + "".join(sam_line(q, p, c, s) for q, p, c, s in specs))
    print(f"wrote {fa} and {sam} ({len(specs)} reads)")


if __name__ == "__main__":
    main()
