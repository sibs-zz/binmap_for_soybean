#!/usr/bin/env python3
"""
Annotate QTL IciMapping qtlout.txt with physical coordinates via bin_physical.tsv.

Example:
  python3 annotate_qtlout_physical.py \\
    --qtlout binmap_test5/qtlout.txt \\
    --physical binmap_test5/bin_physical.tsv \\
    --out binmap_test5/qtlout.with_physical.tsv
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple


def load_physical(path: Path) -> Tuple[Dict[str, dict], Dict[str, List[Tuple[float, int, int, str]]]]:
    phys: Dict[str, dict] = {}
    chrom_bins: Dict[str, List[Tuple[float, int, int, str]]] = {}
    with open(path) as fh:
        header = fh.readline()
        if not header.lower().startswith("bin_id"):
            raise SystemExit(f"Unexpected physical table header in {path}: {header!r}")
        for line in fh:
            if not line.strip():
                continue
            bin_id, chrom, start, end, mid_bp, cm = line.rstrip("\n").split("\t")
            rec = {
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "mid_bp": int(mid_bp),
                "cM": float(cm),
            }
            phys[bin_id] = rec
            chrom_bins.setdefault(chrom, []).append((rec["cM"], rec["start"], rec["end"], bin_id))
    for c in chrom_bins:
        chrom_bins[c].sort()
    return phys, chrom_bins


def ci_physical(
    chrom_bins: Dict[str, List[Tuple[float, int, int, str]]],
    chrom: str,
    left_ci: float,
    right_ci: float,
) -> Tuple[Optional[str], Optional[str], Optional[int], Optional[int]]:
    """Union of bins on chrom whose cM falls within [left_ci, right_ci]."""
    rows = chrom_bins.get(chrom, [])
    if not rows:
        return None, None, None, None
    lo_ci, hi_ci = (left_ci, right_ci) if left_ci <= right_ci else (right_ci, left_ci)
    hits = [r for r in rows if lo_ci - 1e-9 <= r[0] <= hi_ci + 1e-9]
    if not hits:
        lo = min(rows, key=lambda x: abs(x[0] - lo_ci))
        hi = min(rows, key=lambda x: abs(x[0] - hi_ci))
        a, b = (lo, hi) if lo[0] <= hi[0] else (hi, lo)
        hits = [r for r in rows if a[0] - 1e-9 <= r[0] <= b[0] + 1e-9]
    if not hits:
        return None, None, None, None
    return hits[0][3], hits[-1][3], hits[0][1], hits[-1][2]


def parse_qtlout(path: Path) -> List[dict]:
    rows: List[dict] = []
    with open(path) as fh:
        header = fh.readline()
        if not header.strip():
            raise SystemExit(f"Empty qtlout: {path}")
        for line in fh:
            if not line.strip():
                continue
            parts = line.split()
            if len(parts) < 12:
                raise SystemExit(f"Bad qtlout line (need ≥12 fields): {line!r}")
            rows.append(
                {
                    "TraitID": parts[0],
                    "TraitName": parts[1],
                    "Chromosome": parts[2],
                    "Position": parts[3],
                    "LeftMarker": parts[4],
                    "RightMarker": parts[5],
                    "LOD": parts[6],
                    "PVE(%)": parts[7],
                    "Add": parts[8],
                    "Dom": parts[9],
                    "LeftCI": parts[10],
                    "RightCI": parts[11],
                }
            )
    return rows


def annotate(qtlout: Path, physical: Path, out: Path) -> None:
    phys, chrom_bins = load_physical(physical)
    rows = parse_qtlout(qtlout)

    extra_header = [
        "PhysChrom",
        "PeakLeft_start",
        "PeakLeft_end",
        "PeakLeft_mid_bp",
        "PeakRight_start",
        "PeakRight_end",
        "PeakRight_mid_bp",
        "Peak_interval_start",
        "Peak_interval_end",
        "Peak_interval_Mb",
        "CI_LeftBin",
        "CI_RightBin",
        "CI_start",
        "CI_end",
        "CI_Mb",
    ]

    missing = 0
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fh:
        fh.write(
            "\t".join(
                [
                    "TraitID",
                    "TraitName",
                    "Chromosome",
                    "Position_cM",
                    "LeftMarker",
                    "RightMarker",
                    "LOD",
                    "PVE(%)",
                    "Add",
                    "Dom",
                    "LeftCI",
                    "RightCI",
                ]
                + extra_header
            )
            + "\n"
        )
        for r in rows:
            L = phys.get(r["LeftMarker"])
            R = phys.get(r["RightMarker"])
            if not L or not R:
                missing += 1
                vals = [""] * len(extra_header)
            else:
                phys_chrom = L["chrom"]
                peak_start = min(L["start"], R["start"])
                peak_end = max(L["end"], R["end"])
                ci_lb, ci_rb, ci_s, ci_e = ci_physical(
                    chrom_bins, phys_chrom, float(r["LeftCI"]), float(r["RightCI"])
                )
                vals = [
                    phys_chrom,
                    str(L["start"]),
                    str(L["end"]),
                    str(L["mid_bp"]),
                    str(R["start"]),
                    str(R["end"]),
                    str(R["mid_bp"]),
                    str(peak_start),
                    str(peak_end),
                    f"{peak_start / 1e6:.3f}-{peak_end / 1e6:.3f}",
                    ci_lb or "",
                    ci_rb or "",
                    "" if ci_s is None else str(ci_s),
                    "" if ci_e is None else str(ci_e),
                    "" if ci_s is None else f"{ci_s / 1e6:.3f}-{ci_e / 1e6:.3f}",
                ]
            fh.write(
                "\t".join(
                    [
                        r["TraitID"],
                        r["TraitName"],
                        r["Chromosome"],
                        r["Position"],
                        r["LeftMarker"],
                        r["RightMarker"],
                        r["LOD"],
                        r["PVE(%)"],
                        r["Add"],
                        r["Dom"],
                        r["LeftCI"],
                        r["RightCI"],
                    ]
                    + vals
                )
                + "\n"
            )

    print(f"Wrote {out}  (QTL={len(rows)}, missing_bins={missing})")


def main() -> int:
    p = argparse.ArgumentParser(
        description="Add physical coordinates to IciMapping qtlout.txt using bin_physical.tsv"
    )
    p.add_argument("--qtlout", type=Path, required=True, help="IciMapping qtlout.txt")
    p.add_argument("--physical", type=Path, required=True, help="bin_physical.tsv from binmap_pipeline")
    p.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output TSV (default: <qtlout>.with_physical.tsv)",
    )
    args = p.parse_args()
    out = args.out
    if out is None:
        out = args.qtlout.with_suffix(args.qtlout.suffix + ".with_physical.tsv")
        # qtlout.txt -> qtlout.txt.with_physical.tsv is ugly; prefer sibling name
        out = args.qtlout.parent / (args.qtlout.stem + ".with_physical.tsv")
    annotate(args.qtlout, args.physical, out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
