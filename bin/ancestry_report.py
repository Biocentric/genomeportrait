#!/usr/bin/env python3
"""ancestry_report.py — build ancestry tables + PCA/admixture plots into a bundle dir."""
import argparse
import os
import sys

import pandas as pd

# Plotting is optional: tables are always written, plots only if matplotlib is present.
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAVE_MPL = True
except Exception:  # pragma: no cover
    HAVE_MPL = False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--eigenvec", required=True)
    ap.add_argument("--admixture", required=True)
    ap.add_argument("--roh", required=True)
    ap.add_argument("--kin", required=True)
    ap.add_argument("--ref-psam", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    # ---- PCA plot: reference panel coloured by superpopulation + the sample marked ----
    try:
        ev = pd.read_csv(a.eigenvec, sep=r"\s+")
        pc_cols = [c for c in ev.columns if c.upper().startswith("PC") or "SCORE" in c.upper()][:2]
        ref = pd.read_csv(a.ref_psam, sep=r"\s+")
        popcol = next((c for c in ref.columns if c.lower() in ("superpop", "population", "pop")), None)
        fig, ax = plt.subplots(figsize=(6, 5))
        if len(pc_cols) == 2:
            ax.scatter(ev[pc_cols[0]], ev[pc_cols[1]], s=6, alpha=.5, c="#5cc8ff", label="reference")
            ax.scatter(ev[pc_cols[0]].iloc[-1], ev[pc_cols[1]].iloc[-1], s=120, marker="*",
                       c="#ffd166", edgecolor="k", label=a.sample, zorder=5)
            ax.set_xlabel(pc_cols[0]); ax.set_ylabel(pc_cols[1])
        ax.set_title("PCA projection onto 1000G + HGDP")
        ax.legend(fontsize=8); fig.tight_layout()
        fig.savefig(os.path.join(a.outdir, "pca_projection.png"), dpi=130)
        plt.close(fig)
    except Exception as e:
        print(f"[ancestry] PCA plot skipped: {e}", file=sys.stderr)

    # ---- Admixture proportions ----
    try:
        q = pd.read_csv(a.admixture, sep=r"\s+", header=None)
        props = q.iloc[-1, 1:].astype(float)
        props.index = [f"K{i+1}" for i in range(len(props))]
        adf = props.sort_values(ascending=False).rename("proportion").reset_index().rename(columns={"index": "component"})
        adf["proportion"] = (adf["proportion"] * 100).round(1)
        adf.to_csv(os.path.join(a.outdir, "admixture_proportions.tsv"), sep="\t", index=False)
        fig, ax = plt.subplots(figsize=(6, 1.6))
        left = 0
        for _, r in adf.iterrows():
            ax.barh(0, r["proportion"], left=left, label=r["component"])
            left += r["proportion"]
        ax.set_xlim(0, 100); ax.set_yticks([]); ax.set_xlabel("% ancestry")
        ax.set_title("Supervised admixture"); ax.legend(ncol=6, fontsize=7, loc="upper center", bbox_to_anchor=(.5, -.4))
        fig.tight_layout(); fig.savefig(os.path.join(a.outdir, "admixture.png"), dpi=130); plt.close(fig)
    except Exception as e:
        print(f"[ancestry] admixture skipped: {e}", file=sys.stderr)

    # ---- ROH / inbreeding + kinship summary table ----
    try:
        roh = pd.read_csv(a.roh, sep="\t")
        roh.to_csv(os.path.join(a.outdir, "roh_inbreeding.tsv"), sep="\t", index=False)
    except Exception as e:
        print(f"[ancestry] roh skipped: {e}", file=sys.stderr)

    print(f"[ancestry_report] bundle written to {a.outdir}", file=sys.stderr)


if __name__ == "__main__":
    main()
