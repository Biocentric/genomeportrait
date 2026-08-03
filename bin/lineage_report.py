#!/usr/bin/env python3
"""lineage_report.py — combine mtDNA (HaploGrep2) and Y (Yleaf) haplogroup calls."""
import argparse
import sys

import pandas as pd


def mt_haplogroup(path):
    try:
        df = pd.read_csv(path, sep="\t")
        col = next((c for c in df.columns if c.lower() in ("haplogroup", "haplo")), df.columns[1])
        qcol = next((c for c in df.columns if "quality" in c.lower()), None)
        hg = str(df[col].iloc[0])
        q = f" (quality {df[qcol].iloc[0]})" if qcol else ""
        return hg + q
    except Exception as e:
        return f"NA ({e})"


def y_haplogroup(path):
    if not path:
        return "NA (no Y chromosome)"
    try:
        df = pd.read_csv(path, sep="\t")
        col = next((c for c in df.columns if "haplo" in c.lower()), df.columns[-1])
        return str(df[col].iloc[0])
    except Exception as e:
        return f"NA ({e})"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--mtdna", required=True)
    ap.add_argument("--ydna", default=None)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    rows = [
        {"lineage": "Maternal (mtDNA)", "haplogroup": mt_haplogroup(a.mtdna),
         "meaning": "Direct mother-line ancestry (PhyloTree)"},
        {"lineage": "Paternal (Y-DNA)", "haplogroup": y_haplogroup(a.ydna),
         "meaning": "Direct father-line ancestry (ISOGG/Yleaf)"},
    ]
    df = pd.DataFrame(rows)
    with open(a.out, "w") as fh:
        fh.write(f"# lineage for {a.sample}\n")
    df.to_csv(a.out, sep="\t", index=False, mode="a")
    print("[lineage_report] done", file=sys.stderr)


if __name__ == "__main__":
    main()
