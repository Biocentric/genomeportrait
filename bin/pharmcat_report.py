#!/usr/bin/env python3
"""pharmcat_report.py — flatten the PharmCAT JSON report into a gene/phenotype table."""
import argparse
import json
import sys

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--json", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    rows = []
    try:
        with open(a.json) as fh:
            data = json.load(fh)
        genes = data.get("genes", {})
        # PharmCAT groups by source (CPIC / DPWG); iterate defensively
        for source, gmap in genes.items() if isinstance(genes, dict) else []:
            for gene, info in (gmap or {}).items():
                diplotypes = info.get("recommendationDiplotypes") or info.get("sourceDiplotypes") or []
                for dp in diplotypes:
                    rows.append({
                        "source": source,
                        "gene": gene,
                        "diplotype": dp.get("label", dp.get("name", ".")),
                        "phenotype": "; ".join(dp.get("phenotypes", [])) or ".",
                        "activity": dp.get("activityScore", "."),
                    })
    except Exception as e:
        rows.append({"source": "error", "gene": "", "diplotype": str(e), "phenotype": "", "activity": ""})

    df = pd.DataFrame(rows).drop_duplicates()
    with open(a.out, "w") as fh:
        fh.write(f"# pharmacogenomics (pgx) for {a.sample} — informational, not a prescription\n")
    df.to_csv(a.out, sep="\t", index=False, mode="a")
    print(f"[pharmcat_report] {len(df)} gene-phenotype rows", file=sys.stderr)


if __name__ == "__main__":
    main()
