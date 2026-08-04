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
        if len(df) == 0:
            return "NA (no call)"
        # Yleaf's hg_prediction.hg columns are:
        #   Sample_name  Hg  Hg_marker  Total_reads  Valid_markers  QC-score  QC-1..3
        # None of them contain "haplo", so match "Hg" exactly first — otherwise the old
        # fallback to the last column returned the QC score as the haplogroup.
        col = next((c for c in df.columns if c.strip().lower() == "hg"), None)
        if col is None:
            col = next((c for c in df.columns if "haplo" in c.lower()), df.columns[-1])
        hg = str(df[col].iloc[0])
        marker = next((c for c in df.columns if c.strip().lower() == "hg_marker"), None)
        return f"{hg} ({df[marker].iloc[0]})" if marker else hg
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
