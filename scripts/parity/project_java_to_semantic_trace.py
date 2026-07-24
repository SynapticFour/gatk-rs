#!/usr/bin/env python3
"""Project Java HC dump artifacts into semantic-trace NDJSON (v1).

Does not modify GATK algorithms. Builds comparable checkpoints from:
  - Active-region TSV (contig/start/end/is_active[/extended...])
  - Optional assembly-graph / haplotype TSV dumps when present
  - VCF emission rows

Output uses schema gatk_rs.hc.semantic_trace/v1 with impl=java so
compare_semantic_trace.py can locate the first Java↔Rust divergence.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Iterable

SCHEMA = "gatk_rs.hc.semantic_trace/v1"


def parse_bool(s: str) -> bool:
    return s.strip().lower() in {"1", "true", "t", "yes", "y"}


def write_event(
    out: list[dict[str, Any]],
    seq: int,
    stage: str,
    region: dict[str, Any] | None,
    payload: dict[str, Any],
) -> int:
    out.append(
        {
            "schema": SCHEMA,
            "seq": seq,
            "impl": "java",
            "stage": stage,
            "region": region,
            "payload": payload,
        }
    )
    return seq + 1


def load_regions_tsv(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    header: list[str] | None = None
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        if header is None:
            # Accept either header row or bare contig/start/end/is_active
            if cols[0].lower() in {"contig", "chrom", "chromosome"}:
                header = [c.lower() for c in cols]
                continue
            header = ["contig", "start", "end", "is_active"]
        while len(cols) < len(header):
            cols.append("")
        rec = {header[i]: cols[i] for i in range(len(header))}
        contig = rec.get("contig") or rec.get("chrom") or rec.get("chromosome")
        start = int(float(rec["start"]))
        end = int(float(rec["end"]))
        is_active = parse_bool(rec.get("is_active", "true"))
        row = {
            "contig": contig,
            "start": start,
            "end": end,
            "is_active": is_active,
        }
        if "extended_start" in rec and rec["extended_start"]:
            row["extended_start"] = int(float(rec["extended_start"]))
        if "extended_end" in rec and rec["extended_end"]:
            row["extended_end"] = int(float(rec["extended_end"]))
        if "read_count" in rec and rec["read_count"]:
            row["read_count"] = int(float(rec["read_count"]))
        rows.append(row)
    return rows


def load_vcf_sites(path: Path) -> list[dict[str, Any]]:
    sites: list[dict[str, Any]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 5:
            continue
        chrom, pos, _id, ref, alt = cols[:5]
        qual_s = cols[5] if len(cols) > 5 else "."
        filt = cols[6] if len(cols) > 6 else "."
        qual = None if qual_s in {".", ""} else float(qual_s)
        sites.append(
            {
                "chrom": chrom,
                "pos": int(pos),
                "ref": ref,
                "alt": alt.split(",") if alt != "." else [],
                "qual": qual,
                "filter": [] if filt in {".", "PASS"} else filt.split(";"),
            }
        )
    sites.sort(key=lambda s: (s["chrom"], s["pos"], s["ref"], ",".join(s["alt"])))
    return sites


def assign_sites_to_regions(
    sites: list[dict[str, Any]], regions: list[dict[str, Any]]
) -> dict[tuple[Any, ...], list[dict[str, Any]]]:
    by_region: dict[tuple[Any, ...], list[dict[str, Any]]] = {
        (r["contig"], r["start"], r["end"]): [] for r in regions
    }
    unmatched: list[dict[str, Any]] = []
    for site in sites:
        placed = False
        for r in regions:
            if site["chrom"] != r["contig"]:
                continue
            if r["start"] <= site["pos"] <= r["end"]:
                by_region[(r["contig"], r["start"], r["end"])].append(site)
                placed = True
                break
        if not placed:
            unmatched.append(site)
    if unmatched:
        by_region[("__unmatched__", 0, 0)] = unmatched
    return by_region


def project(
    regions: list[dict[str, Any]],
    sites: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    events: list[dict[str, Any]] = []
    seq = 1
    for r in regions:
        region = {
            "contig": r["contig"],
            "start": r["start"],
            "end": r["end"],
            "is_active": r["is_active"],
        }
        # Activity profile cut (Java dump often equals active-region span).
        seq = write_event(
            events,
            seq,
            "activity_profile",
            region,
            {
                "padded_start": r.get("extended_start", r["start"]),
                "padded_end": r.get("extended_end", r["end"]),
                "extension": 0,
                "span_bp": r["end"] - r["start"] + 1,
            },
        )
        seq = write_event(
            events,
            seq,
            "active_region",
            region,
            {
                "is_active": r["is_active"],
                "extended_start": r.get("extended_start", r["start"]),
                "extended_end": r.get("extended_end", r["end"]),
                "read_count": r.get("read_count", 0),
                "pileup_locus_count": 0,
            },
        )

    by_region = assign_sites_to_regions(sites, regions)
    for r in regions:
        key = (r["contig"], r["start"], r["end"])
        region_sites = by_region.get(key, [])
        if not r["is_active"]:
            seq = write_event(
                events,
                seq,
                "inactive_rcm",
                {
                    "contig": r["contig"],
                    "start": r["start"],
                    "end": r["end"],
                    "is_active": False,
                },
                {"locus_count": 0},
            )
        region = {
            "contig": r["contig"],
            "start": r["start"],
            "end": r["end"],
            "is_active": r["is_active"],
        }
        seq = write_event(
            events,
            seq,
            "vcf_emission",
            region,
            {"record_count": len(region_sites), "sites": region_sites},
        )

    unmatched = by_region.get(("__unmatched__", 0, 0))
    if unmatched:
        seq = write_event(
            events,
            seq,
            "vcf_emission",
            None,
            {"record_count": len(unmatched), "sites": unmatched},
        )
    return events


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--regions",
        type=Path,
        help="Java/Rust assembly-regions TSV (contig start end is_active ...)",
    )
    p.add_argument("--vcf", type=Path, help="Java HaplotypeCaller VCF")
    p.add_argument("-o", "--output", required=True, type=Path)
    args = p.parse_args()

    regions = load_regions_tsv(args.regions) if args.regions else []
    sites = load_vcf_sites(args.vcf) if args.vcf else []
    if not regions and not sites:
        raise SystemExit("provide --regions and/or --vcf")

    # If only VCF: synthesize one synthetic active region per contig span.
    if not regions and sites:
        by_chrom: dict[str, list[dict[str, Any]]] = {}
        for s in sites:
            by_chrom.setdefault(s["chrom"], []).append(s)
        for chrom, ss in sorted(by_chrom.items()):
            lo = min(x["pos"] for x in ss)
            hi = max(x["pos"] for x in ss)
            regions.append(
                {
                    "contig": chrom,
                    "start": lo,
                    "end": hi,
                    "is_active": True,
                }
            )

    events = project(regions, sites)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8") as fh:
        for ev in events:
            fh.write(json.dumps(ev, separators=(",", ":"), sort_keys=True))
            fh.write("\n")
    print(f"wrote {len(events)} events → {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
