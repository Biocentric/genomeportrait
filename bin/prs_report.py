#!/usr/bin/env python3
"""prs_report.py — turn PLINK2 --score sums into per-trait polygenic scores."""
import argparse
import sys

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--sscore", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--min-overlap", type=float, default=0.75)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    # PLINK2 .sscore: per-PGS columns named <PGS>_SUM ; ALLELE_CT / NAMED_ALLELE_DOSAGE_SUM track overlap
    sc = pd.read_csv(a.sscore, sep=r"\s+")
    try:
        meta = pd.read_csv(a.metadata, sep="\t")
        labels = dict(zip(meta["pgs_id"], meta["trait"]))
    except Exception:
        labels = {}

    rows = []
    score_cols = [c for c in sc.columns if c.endswith("_SUM") or "SCORE" in c.upper()]
    for c in score_cols:
        pgs = c.replace("_SUM", "").replace("_AVG", "")
        val = float(sc[c].iloc[0]) if len(sc) else float("nan")
        rows.append({
            "pgs_id": pgs,
            "trait": labels.get(pgs, "see PGS Catalog"),
            "raw_score": round(val, 4),
            "note": "additive sum of effect-allele dosages × weights",
        })

    df = pd.DataFrame(rows)
    if len(df):
        # Z-standardise across the requested scores so the report can rank them
        m, s = df["raw_score"].mean(), df["raw_score"].std(ddof=0) or 1.0
        df["z_relative"] = ((df["raw_score"] - m) / s).round(2)
    with open(a.out, "w") as fh:
        fh.write(f"# polygenic prs for {a.sample} (relative ranking only; not absolute risk)\n")
    df.to_csv(a.out, sep="\t", index=False, mode="a")
    print(f"[prs_report] {len(df)} scores written", file=sys.stderr)


if __name__ == "__main__":
    main()
