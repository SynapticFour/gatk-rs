#!/usr/bin/env python3
"""Extract L4 P12 Java FORMAT fixtures from the pinned Phase E Java VCF.

Source (default): parity/reports/p12_realworld_na12878_20k.java.vcf
Output: parity/fixtures/p12-java-format/all_sites.tsv (+ optional per-site TSVs)

Regenerate after refreshing the Java baseline:
  ./scripts/parity/run_p12_realworld_na12878_20k.sh   # or your Java HC refresh
  python3 scripts/parity/extract_p12_java_format_fixtures.py
"""

from __future__ import annotations

import argparse
import pathlib
import sys


def parse_info(info: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for part in info.split(";"):
        if "=" not in part:
            continue
        k, v = part.split("=", 1)
        out[k] = v
    return out


def parse_format_sample(format_keys: str, sample: str) -> dict[str, str]:
    keys = format_keys.split(":")
    vals = sample.split(":")
    if len(keys) != len(vals):
        raise ValueError(f"FORMAT/SAMPLE length mismatch: {format_keys!r} {sample!r}")
    return dict(zip(keys, vals))


def extract(vcf_path: pathlib.Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with vcf_path.open() as fh:
        for line in fh:
            if not line.startswith("2\t"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 10:
                continue
            chrom, pos, _id, ref, alt, qual, _filter, info, fmt_keys, sample = cols[:10]
            fmt = parse_format_sample(fmt_keys, sample)
            info_d = parse_info(info)
            af = info_d.get("AF", ".")
            if "," in af:
                af = af.split(",")[0]
            rows.append(
                {
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "gt": fmt.get("GT", "."),
                    "pl": fmt.get("PL", "."),
                    "gq": fmt.get("GQ", "."),
                    "ad": fmt.get("AD", "."),
                    "dp": fmt.get("DP", "."),
                    "qual": qual,
                    "af": af,
                }
            )
    rows.sort(key=lambda r: int(r["pos"]))
    return rows


def write_all_sites(out_dir: pathlib.Path, rows: list[dict[str, str]]) -> pathlib.Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "all_sites.tsv"
    header = ["chrom", "pos", "ref", "alt", "gt", "pl", "gq", "ad", "dp", "qual", "af"]
    with path.open("w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(r[k] for k in header) + "\n")
    return path


def write_per_site(out_dir: pathlib.Path, rows: list[dict[str, str]]) -> int:
    sites_dir = out_dir / "sites"
    sites_dir.mkdir(parents=True, exist_ok=True)
    header = ["pos", "ref", "alt", "gt", "pl", "gq", "ad", "dp", "qual", "af"]
    for r in rows:
        name = f"{r['pos']}_{r['ref']}_{r['alt']}.tsv"
        path = sites_dir / name
        with path.open("w") as fh:
            fh.write("\t".join(header) + "\n")
            fh.write("\t".join(r[k] for k in header) + "\n")
    return len(rows)


def main() -> int:
    repo = pathlib.Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--vcf",
        type=pathlib.Path,
        default=repo / "parity/reports/p12_realworld_na12878_20k.java.vcf",
        help="Pinned Java P12 VCF",
    )
    parser.add_argument(
        "--out-dir",
        type=pathlib.Path,
        default=repo / "parity/fixtures/p12-java-format",
    )
    parser.add_argument(
        "--per-site",
        action="store_true",
        help="Also write parity/fixtures/p12-java-format/sites/<pos>_<ref>_<alt>.tsv",
    )
    args = parser.parse_args()
    if not args.vcf.is_file():
        print(f"error: VCF not found: {args.vcf}", file=sys.stderr)
        return 1
    rows = extract(args.vcf)
    if len(rows) != 66:
        print(
            f"error: expected 66 variant rows, got {len(rows)} from {args.vcf}",
            file=sys.stderr,
        )
        return 1
    all_path = write_all_sites(args.out_dir, rows)
    print(f"wrote {all_path} ({len(rows)} rows)")
    if args.per_site:
        n = write_per_site(args.out_dir, rows)
        print(f"wrote {n} per-site TSVs under {args.out_dir / 'sites'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
