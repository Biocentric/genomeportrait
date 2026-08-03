#!/usr/bin/env python3
"""extras_report.py — combine telomere-length estimate and mtDNA heteroplasmy count."""
import argparse
import sys

import pandas as pd


def telomere_len(path):
    if not path:
        return None
    try:
        df = pd.read_csv(path, sep="\t", comment="#")
        col = next((c for c in df.columns if "LENGTH" in c.upper() or "tel" in c.lower()), None)
        if col:
            return round(float(df[col].dropna().iloc[0]), 2)
    except Exception as e:
        print(f"[extras] telseq parse: {e}", file=sys.stderr)
    return None


def heteroplasmy_count(path):
    if not path:
        return None
    try:
        df = pd.read_csv(path, sep="\t")
        lvl = next((c for c in df.columns if "level" in c.lower()), None)
        if lvl:
            sub = df[(df[lvl].astype(float) >= 0.01) & (df[lvl].astype(float) <= 0.95)]
            return int(len(sub))
        return int(len(df))
    except Exception as e:
        print(f"[extras] mutserve parse: {e}", file=sys.stderr)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--telseq", default=None)
    ap.add_argument("--mito", default=None)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    rows = []
    tl = telomere_len(a.telseq)
    if tl is not None:
        rows.append({"metric": "Mean telomere length", "value": f"{tl} kb",
                     "note": "estimated from WGS read content (Telseq)"})
    hc = heteroplasmy_count(a.mito)
    if hc is not None:
        rows.append({"metric": "mtDNA heteroplasmies (1–95%)", "value": hc,
                     "note": "low-level mitochondrial variants (mutserve)"})
    if not rows:
        rows.append({"metric": "extras", "value": "none", "note": "modules skipped or no data"})

    df = pd.DataFrame(rows)
    with open(a.out, "w") as fh:
        fh.write(f"# extras for {a.sample}\n")
    df.to_csv(a.out, sep="\t", index=False, mode="a")
    print("[extras_report] done", file=sys.stderr)


if __name__ == "__main__":
    main()
