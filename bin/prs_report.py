#!/usr/bin/env python3
"""prs_report.py — collate per-PGS PLINK2 --score outputs into a polygenic-score table."""
import argparse
import bisect
import glob
import os
import sys

import pandas as pd



def load_reference(dirpath):
    """Reference-panel score distributions: {pgs_id: {group: sorted list of scores}}.

    A polygenic score only means something relative to a population scored with the SAME
    file, and PGS distributions differ strongly by ancestry, so keep them per superpopulation
    rather than pooling."""
    dist = {}
    if not dirpath or not os.path.isdir(dirpath):
        return dist
    psam = os.path.join(dirpath, "reference.psam")
    groups = {}
    if os.path.exists(psam):
        try:
            ref = pd.read_csv(psam, sep=r"\s+")
            ref.columns = [c.lstrip("#") for c in ref.columns]
            idc = ref.columns[0]
            gcol = next((c for c in ref.columns if c.lower() in ("superpop", "population", "pop")), None)
            if gcol:
                groups = dict(zip(ref[idc].astype(str), ref[gcol].astype(str)))
        except Exception as e:
            print(f"[prs_report] reference psam: {e}", file=sys.stderr)
    for f in sorted(glob.glob(os.path.join(dirpath, "*.sscore"))):
        pgs_id = os.path.basename(f)[:-7].split("_")[0]
        try:
            d = pd.read_csv(f, sep=r"\s+")
        except Exception:
            continue
        sum_col = next((c for c in d.columns if c.upper().endswith("_SUM") and "DOSAGE" not in c.upper()), None)
        idc = next((c for c in d.columns if c.upper().lstrip("#") in ("IID", "FID")), d.columns[0])
        if not sum_col:
            continue
        by = {"ALL": []}
        for iid, val in zip(d[idc].astype(str), pd.to_numeric(d[sum_col], errors="coerce")):
            if pd.isna(val):
                continue
            by["ALL"].append(float(val))
            g = groups.get(iid)
            if g:
                by.setdefault(g, []).append(float(val))
        dist[pgs_id] = {g: sorted(v) for g, v in by.items() if len(v) >= 20}
    return dist


def percentile_of(value, sorted_vals):
    if not sorted_vals or value is None or pd.isna(value):
        return None
    below = bisect.bisect_left(sorted_vals, float(value))
    return round(100.0 * below / len(sorted_vals), 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--indir", default=".")
    ap.add_argument("--metadata", default=None)
    ap.add_argument("--ref-scores", default=None,
                    help="dir of reference-panel .sscore files + reference.psam")
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

    refdist = load_reference(a.ref_scores)
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
        # plink2's DENOM/ALLELE_CT is an ALLELE count (2 per diploid variant), not a variant
        # count — dividing by 2 gives the variants actually used.
        alleles = int(df[nvar_col].iloc[0]) if nvar_col else None
        used = alleles // 2 if alleles is not None else None
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
            "alleles_counted": alleles if alleles is not None else "NA",
            "variants_in_score": total,
            "coverage": cover,
            "pgs_catalog": ex.get("pgs_catalog", f"https://www.pgscatalog.org/score/{pgs_id}/"),
        })
        raw = rows[-1]["raw_score"]
        for grp, vals in sorted((refdist.get(pgs_id) or {}).items()):
            pct = percentile_of(raw, vals)
            if pct is not None:
                rows[-1][f"percentile_{grp}"] = pct

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
