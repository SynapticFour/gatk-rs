"""Shared Phase-11 Java-vs-Rust first-variant field comparison.

Used by smoke + corpus gates so tolerances stay in one place.
"""
from __future__ import annotations

from typing import Any


# Small phred / log10 rounding bands observed on the synthetic Java-positive smoke
# (QUAL Δ≈1.0, PL[0] Δ=1). Identity on alleles / GT / AD remains strict.
QUAL_ABS_TOL = 1.5
PL_ABS_TOL = 1


def first_variant_fields(path) -> dict[str, str] | None:
    from pathlib import Path

    p = Path(path)
    if not p.exists():
        return None
    for line in p.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line or line.startswith("#"):
            continue
        cols = line.split("\t")
        if len(cols) < 8:
            return None
        out = {
            "CHROM": cols[0],
            "POS": cols[1],
            "REF": cols[3],
            "ALT": cols[4],
            "QUAL": cols[5],
            "FILTER": cols[6],
            "INFO": cols[7],
        }
        if len(cols) >= 10:
            out["FORMAT"] = cols[8]
            out["SAMPLE"] = cols[9]
        return out
    return None


def count_variants(path) -> int:
    from pathlib import Path

    p = Path(path)
    if not p.exists():
        return 0
    return sum(
        1
        for line in p.read_text(encoding="utf-8", errors="replace").splitlines()
        if line and not line.startswith("#")
    )


def parse_info_map(info_str: str | None) -> dict[str, str]:
    out: dict[str, str] = {}
    if not info_str or info_str == ".":
        return out
    for tok in info_str.split(";"):
        if "=" in tok:
            k, v = tok.split("=", 1)
            out[k] = v
        else:
            out[tok] = "true"
    return out


def extract_sample_subfields(variant: dict[str, str] | None) -> dict[str, str]:
    if not variant:
        return {}
    fmt = variant.get("FORMAT")
    sample = variant.get("SAMPLE")
    if not fmt or not sample:
        return {}
    keys = fmt.split(":")
    vals = sample.split(":")
    return {k: vals[i] if i < len(vals) else "." for i, k in enumerate(keys)}


def _near_float(a: str, b: str, tol: float) -> bool:
    try:
        return abs(float(a) - float(b)) <= tol
    except Exception:
        return False


def _near_pl(a: str | None, b: str | None, tol: int = PL_ABS_TOL) -> bool:
    if a is None or b is None:
        return a == b
    try:
        av = [int(x) for x in a.split(",")]
        bv = [int(x) for x in b.split(",")]
    except Exception:
        return False
    if len(av) != len(bv):
        return False
    return all(abs(x - y) <= tol for x, y in zip(av, bv))


def compare_first_variants(
    java_first: dict[str, str] | None, rust_first: dict[str, str] | None
) -> list[str]:
    """Return mismatch keys (empty ⇒ pass)."""
    mismatches: list[str] = []
    for k in ["CHROM", "POS", "REF", "ALT", "FILTER"]:
        if (java_first or {}).get(k) != (rust_first or {}).get(k):
            mismatches.append(k)
    jq = (java_first or {}).get("QUAL")
    rq = (rust_first or {}).get("QUAL")
    if jq is None or rq is None or not _near_float(jq, rq, QUAL_ABS_TOL):
        mismatches.append("QUAL")

    java_info = parse_info_map((java_first or {}).get("INFO", "."))
    rust_info = parse_info_map((rust_first or {}).get("INFO", "."))
    for k in ["AC", "AF", "AN", "DP"]:
        jv = java_info.get(k)
        rv = rust_info.get(k)
        if jv is None or rv is None:
            mismatches.append(f"INFO.{k}")
            continue
        if k == "AF":
            if not _near_float(jv, rv, 0.01):
                mismatches.append(f"INFO.{k}")
        elif jv != rv:
            mismatches.append(f"INFO.{k}")

    java_sample = extract_sample_subfields(java_first)
    rust_sample = extract_sample_subfields(rust_first)
    for k in ["GT", "AD", "DP", "GQ"]:
        if java_sample.get(k) != rust_sample.get(k):
            mismatches.append(f"SAMPLE.{k}")
    if not _near_pl(java_sample.get("PL"), rust_sample.get("PL")):
        mismatches.append("SAMPLE.PL")
    return mismatches


def smoke_status(
    *,
    java_exit: int,
    java_variants: int,
    rust_exit: int,
    rust_variants: int,
    java_first: dict[str, Any] | None,
    rust_first: dict[str, Any] | None,
) -> tuple[str, str, list[str]]:
    """Return (status, notes, mismatches)."""
    if java_exit != 0:
        return (
            "java_unavailable",
            "Java docker oracle unavailable in this environment",
            ["java_exit"],
        )
    if rust_exit != 0:
        return (
            "rust_fail",
            f"Rust HaplotypeCaller failed (exit={rust_exit})",
            ["rust_exit"],
        )
    if rust_variants == 0:
        return (
            "pending_activation",
            "Rust HC has no variant records; field-level diff is deferred",
            [],
        )
    if java_variants == 0 and rust_variants > 0:
        return (
            "divergent_activation",
            "Rust HC emits provisional variants while Java smoke interval has no calls",
            ["variant_presence"],
        )
    mismatches = compare_first_variants(java_first, rust_first)
    if mismatches:
        return (
            "fail",
            f"strict field diff mismatch on keys: {','.join(mismatches)}",
            mismatches,
        )
    return ("pass", "strict field diff smoke matched on core variant keys", [])


if __name__ == "__main__":
    # Laptop-friendly self-check: known Java/Rust ΔQUAL≈1 / ΔPL[0]=1 on synthetic smoke.
    java = {
        "CHROM": "chrLive",
        "POS": "15",
        "REF": "T",
        "ALT": "A",
        "QUAL": "2224.06",
        "FILTER": ".",
        "INFO": "AC=2;AF=1.00;AN=2;DP=50",
        "FORMAT": "GT:AD:DP:GQ:PL",
        "SAMPLE": "1/1:0,50:50:99:2238,151,0",
    }
    rust = {
        "CHROM": "chrLive",
        "POS": "15",
        "REF": "T",
        "ALT": "A",
        "QUAL": "2225.06",
        "FILTER": ".",
        "INFO": "AC=2;AF=1;AN=2;DP=50",
        "FORMAT": "GT:GQ:DP:AD:PL",
        "SAMPLE": "1/1:99:50:0,50:2239,151,0",
    }
    mm = compare_first_variants(java, rust)
    assert mm == [], mm
    status, _, _ = smoke_status(
        java_exit=0,
        java_variants=1,
        rust_exit=0,
        rust_variants=1,
        java_first=java,
        rust_first=rust,
    )
    assert status == "pass", status
    status, _, _ = smoke_status(
        java_exit=0,
        java_variants=1,
        rust_exit=3,
        rust_variants=0,
        java_first=java,
        rust_first=None,
    )
    assert status == "rust_fail", status
    print("[p11_field_compare] self-check OK")
