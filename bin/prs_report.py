#!/usr/bin/env python3
"""prs_report.py — collate per-PGS PLINK2 --score outputs into a polygenic-score table."""
import argparse
import glob
import os
import sys

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--indir", default=".")
    ap.add_argument("--metadata", default=None)
    ap.add_argument("--min-overlap", type=float, default=0.0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    labels = {}
    if a.metadata and os.path.exists(a.metadata):
        try:
            m = pd.read_csv(a.metadata, sep="\t")
            labels = {str(r["pgs_id"]): str(r["trait"]) for _, r in m.iterrows()}
        except Exception:
            pass

    rows = []
    for f in sorted(glob.glob(os.path.join(a.indir, "*.sscore"))):
        acc = os.path.basename(f)[:-7]                      # strip .sscore
        if acc.startswith(a.sample + "."):
            acc = acc[len(a.sample) + 1:]
        pgs_id = acc.split("_")[0]                          # PGS000297_hmPOS_... -> PGS000297
        try:
            df = pd.read_csv(f, sep=r"\s+")
        except Exception:
            continue
        if len(df) == 0:
            continue
        # plink2 names the score columns after the weight column header (BETA_SUM/BETA_AVG here,
        # SCORE1_* by default). Match the *_SUM / *_AVG columns, ignoring the dosage-sum column.
        sum_col = next((c for c in df.columns if c.upper().endswith("_SUM") and "DOSAGE" not in c.upper()), None)
        avg_col = next((c for c in df.columns if c.upper().endswith("_AVG")), None)
        nvar_col = next((c for c in df.columns if c.upper() == "DENOM" or "ALLELE_CT" in c.upper()), None)
        rows.append({
            "pgs_id": pgs_id,
            "trait": labels.get(pgs_id, "see PGS Catalog"),
            "raw_score": round(float(df[sum_col].iloc[0]), 5) if sum_col else float("nan"),
            "per_allele_avg": round(float(df[avg_col].iloc[0]), 8) if avg_col else float("nan"),
            "variants_used": int(df[nvar_col].iloc[0]) if nvar_col else "NA",
        })

    df = pd.DataFrame(rows)
    if len(df) and df["raw_score"].notna().any():
        m, s = df["raw_score"].mean(), df["raw_score"].std(ddof=0) or 1.0
        df["z_relative"] = ((df["raw_score"] - m) / s).round(2)
    with open(a.out, "w") as fh:
        fh.write(f"# polygenic scores for {a.sample} (relative ranking only; not absolute risk)\n")
    (df if len(df) else pd.DataFrame([{"pgs_id": "none", "trait": "no scores computed"}])
     ).to_csv(a.out, sep="\t", index=False, mode="a")
    print(f"[prs_report] {len(df)} scores", file=sys.stderr)


if __name__ == "__main__":
    main()
