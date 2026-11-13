#!/usr/bin/env python3
"""ancestry_report.py — ancestry tables + PCA/admixture plots into a bundle dir.

Inputs come from the BAM_ANCESTRY subworkflow:
  --eigenvec   IID PC1..PC10   (reference + query, projected onto the panel basis)
  --admixture  IID label Q1..QK (supervised ADMIXTURE; label '-' = query, else superpop/pop)
The query is located by IID (== --sample); plink/admixture reorder rows, so never by position.
"""
import argparse
import os
import sys

import pandas as pd

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAVE_MPL = True
except Exception:  # pragma: no cover
    HAVE_MPL = False



def write_proximity(ev, ref, idcol, sample, outdir, n_top=8):
    """Which reference populations does this genome actually sit closest to?

    Distance from the sample to each population's centroid in PC1-4 space (the PCs that carry
    the continental/sub-continental signal). This is the readable version of the PCA."""
    try:
        pcs = [c for c in ("PC1", "PC2", "PC3", "PC4") if c in ev.columns]
        q = ev[ev["IID"] == sample]
        if not pcs or len(q) == 0:
            return
        qv = q.iloc[0][pcs].astype(float).values
        rows = []
        for col in ("Population", "SuperPop"):
            if col not in ref.columns:
                continue
            lab = dict(zip(ref[idcol].astype(str), ref[col].astype(str)))
            d = ev[ev["IID"] != sample].copy()
            d["_g"] = d["IID"].map(lab)
            d = d.dropna(subset=["_g"])
            if len(d) == 0:
                continue
            cent = d.groupby("_g")[pcs].mean()
            dist = ((cent - qv) ** 2).sum(axis=1) ** 0.5
            for grp, dv in dist.sort_values().head(n_top).items():
                rows.append({"level": "population" if col == "Population" else "superpopulation",
                             "closest_group": grp, "distance": round(float(dv), 4),
                             "n_reference_samples": int((d["_g"] == grp).sum())})
        if rows:
            pd.DataFrame(rows).to_csv(os.path.join(outdir, "pca_closest_populations.tsv"),
                                      sep="\t", index=False)
    except Exception as e:
        print(f"[ancestry] proximity summary skipped: {e}", file=sys.stderr)


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

    # reference labels: #IID -> SuperPop / Population
    ref = pd.read_csv(a.ref_psam, sep=r"\s+")
    ref.columns = [c.lstrip("#") for c in ref.columns]
    idcol = ref.columns[0]
    popcol = next((c for c in ref.columns if c.lower() in ("superpop", "population", "pop")), None)
    labels = dict(zip(ref[idcol].astype(str), ref[popcol].astype(str))) if popcol else {}

    # ---- PCA: reference coloured by superpopulation + the query star ----
    try:
        ev = pd.read_csv(a.eigenvec, sep="\t")
        ev["IID"] = ev["IID"].astype(str)
        ev["grp"] = ev["IID"].map(lambda i: labels.get(i, "QUERY" if i == a.sample else "ref"))
        q = ev[ev["IID"] == a.sample]
        # Full coordinates are kept as data but NOT rendered: 3200+ reference rows x 10 PCs is
        # unreadable in a report. The plot carries that; the proximity table is the takeaway.
        ev.to_csv(os.path.join(a.outdir, "pca_coordinates.raw.tsv"), sep="\t", index=False)
        write_proximity(ev, ref, idcol, a.sample, a.outdir)
        if HAVE_MPL and {"PC1", "PC2"}.issubset(ev.columns):
            fig, ax = plt.subplots(figsize=(6, 5))
            for grp, sub in ev[ev["grp"] != "QUERY"].groupby("grp"):
                ax.scatter(sub["PC1"], sub["PC2"], s=7, alpha=.55, label=grp)
            if len(q):
                ax.scatter(q["PC1"], q["PC2"], s=240, marker="*", c="#111", edgecolor="w",
                           linewidth=1.2, label=a.sample, zorder=6)
            ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
            ax.set_title("PCA projection onto 1000G + HGDP")
            ax.legend(fontsize=8, markerscale=1.5)
            fig.tight_layout(); fig.savefig(os.path.join(a.outdir, "pca_projection.png"), dpi=130)
            plt.close(fig)
    except Exception as e:
        print(f"[ancestry] PCA skipped: {e}", file=sys.stderr)

    # ---- Supervised ADMIXTURE: map Q columns -> populations, extract the query ----
    try:
        q = pd.read_csv(a.admixture, sep="\t")
        q["IID"] = q["IID"].astype(str)
        qcols = [c for c in q.columns if c.upper().startswith("Q")]
        lab = q[(q["label"] != "-") & (q["label"].astype(str) != "QUERY")]
        # each population's column = the Q column with the highest mean among its labelled rows
        col2pop = {}
        for pop, sub in lab.groupby("label"):
            col2pop[sub[qcols].mean().idxmax()] = str(pop)
        names = [col2pop.get(c, c) for c in qcols]
        qrow = q[q["IID"] == a.sample]
        if len(qrow):
            frac = qrow.iloc[0][qcols].astype(float).values
            adf = (pd.DataFrame({"population": names, "fraction_pct": (frac * 100).round(2)})
                   .groupby("population", as_index=False)["fraction_pct"].sum()
                   .sort_values("fraction_pct", ascending=False).reset_index(drop=True))
            adf.to_csv(os.path.join(a.outdir, "admixture_proportions.tsv"), sep="\t", index=False)
            if HAVE_MPL:
                fig, ax = plt.subplots(figsize=(6, 1.7))
                left = 0
                for _, r in adf.iterrows():
                    ax.barh(0, r["fraction_pct"], left=left, label=f"{r['population']} {r['fraction_pct']:.0f}%")
                    left += r["fraction_pct"]
                ax.set_xlim(0, 100); ax.set_yticks([]); ax.set_xlabel("% ancestry")
                ax.set_title(f"Supervised admixture — {a.sample}")
                ax.legend(ncol=5, fontsize=7, loc="upper center", bbox_to_anchor=(.5, -.5))
                fig.tight_layout(); fig.savefig(os.path.join(a.outdir, "admixture.png"), dpi=130)
                plt.close(fig)
            top = adf.iloc[0]
            print(f"[ancestry] {a.sample}: {top['population']} {top['fraction_pct']:.1f}%", file=sys.stderr)
    except Exception as e:
        print(f"[ancestry] admixture skipped: {e}", file=sys.stderr)

    # ---- ROH / inbreeding + kinship ----
    for src, name in ((a.roh, "roh_inbreeding.tsv"), (a.kin, "kinship.tsv")):
        try:
            pd.read_csv(src, sep="\t").to_csv(os.path.join(a.outdir, name), sep="\t", index=False)
        except Exception as e:
            print(f"[ancestry] {name} skipped: {e}", file=sys.stderr)

    print(f"[ancestry_report] bundle written to {a.outdir}", file=sys.stderr)


if __name__ == "__main__":
    main()
