#!/usr/bin/env python3
"""pharmcat_report.py — flatten the PharmCAT JSON into a gene/phenotype table.

Accepts either the reporter JSON (top-level "genes") or the phenotyper JSON
(top-level "geneReports"); both nest as {source: {gene: info}}.
"""
import argparse
import json
import sys

import pandas as pd

COLUMNS = ["source", "gene", "diplotype", "phenotype", "activity"]


def diplotype_label(dp):
    """PharmCAT diplotypes carry allele1/allele2 objects rather than a flat label."""
    for key in ("label", "name"):
        if dp.get(key):
            return str(dp[key])
    a1 = (dp.get("allele1") or {}).get("name")
    a2 = (dp.get("allele2") or {}).get("name")
    if a1 and a2:
        return f"{a1}/{a2}"
    return a1 or a2 or "."


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--json", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    rows, note = [], None
    try:
        with open(a.json) as fh:
            data = json.load(fh)
        # reporter JSON -> "genes"; phenotyper JSON -> "geneReports"
        genes = data.get("genes") or data.get("geneReports") or {}
        if not isinstance(genes, dict) or not genes:
            note = "PharmCAT produced no gene reports"
        for source, gmap in (genes.items() if isinstance(genes, dict) else []):
            for gene, info in (gmap or {}).items():
                if not isinstance(info, dict):
                    continue
                diplotypes = info.get("recommendationDiplotypes") or info.get("sourceDiplotypes") or []
                for dp in diplotypes:
                    if not isinstance(dp, dict):
                        continue
                    rows.append({
                        "source": source,
                        "gene": gene or dp.get("gene", "."),
                        "diplotype": diplotype_label(dp),
                        "phenotype": "; ".join(dp.get("phenotypes") or []) or ".",
                        "activity": dp.get("activityScore", "."),
                    })
    except Exception as e:
        note = f"could not parse PharmCAT JSON: {e}"

    df = pd.DataFrame(rows, columns=COLUMNS).drop_duplicates()
    if df.empty:
        # Always emit the header row: an empty frame writes a headerless file, and the final
        # report then fails with pandas' "No columns to parse from file".
        df = pd.DataFrame([{ "source": ".", "gene": ".", "diplotype": ".",
                             "phenotype": note or "no results", "activity": "." }],
                          columns=COLUMNS)
    with open(a.out, "w") as fh:
        fh.write(f"# pharmacogenomics (pgx) for {a.sample} — informational, not a prescription\n")
    df.to_csv(a.out, sep="\t", index=False, mode="a")
    print(f"[pharmcat_report] {len(df)} gene-phenotype rows"
          + (f" ({note})" if note else ""), file=sys.stderr)


if __name__ == "__main__":
    main()
