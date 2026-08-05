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
        ex = meta_extra.get(pgs_id, {})
        used = int(df[nvar_col].iloc[0]) if nvar_col else None
        total = ex.get("n_variants_in_score", ".")
        try:
            cover = f"{100.0 * used / float(total):.1f}%" if used and str(total).isdigit() else "."
        except (ValueError, ZeroDivisionError):
            cover = "."
        rows.append({
            "trait": labels.get(pgs_id, "(trait not stated)"),
            "pgs_id": pgs_id,
            "raw_score": round(float(df[sum_col].iloc[0]), 5) if sum_col else float("nan"),
            "per_allele_avg": round(float(df[avg_col].iloc[0]), 8) if avg_col else float("nan"),
            "variants_used": used if used is not None else "NA",
            "variants_in_score": total,
            "coverage": cover,
            "pgs_catalog": ex.get("pgs_catalog", f"https://www.pgscatalog.org/score/{pgs_id}/"),
        })

    df = pd.DataFrame(rows)
    # NOTE: no cross-score z-score here. The previous version standardised a person's scores
    # against EACH OTHER, which is meaningless — a height score and a BMI score are different
    # units over different variant counts, so their spread says nothing. A PGS is only
    # interpretable as a percentile against a reference population scored with the SAME
    # scoring file, which this pipeline does not yet compute.
    with open(a.out, "w") as fh:
        fh.write(f"# polygenic scores for {a.sample}\n")
        fh.write("# raw_score is the summed weighted dosage; its magnitude depends on how many\n")
        fh.write("# variants the score contains, so scores are NOT comparable with each other\n")
        fh.write("# and a single value carries no meaning without a reference distribution.\n")
    (df if len(df) else pd.DataFrame([{"pgs_id": "none", "trait": "no scores computed"}])
     ).to_csv(a.out, sep="\t", index=False, mode="a")
    print(f"[prs_report] {len(df)} scores", file=sys.stderr)


if __name__ == "__main__":
    main()
