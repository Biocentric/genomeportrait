#!/usr/bin/env python3
"""hla_report.py — tidy T1K genotype output into a per-gene HLA allele table."""
import argparse
import sys

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--t1k", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    # T1K _genotype.tsv columns: gene, num_alleles, allele1, abundance1, quality1, allele2, ...
    rows = []
    try:
        with open(a.t1k) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 3:
                    continue
                gene = f[0]
                a1 = f[2] if len(f) > 2 else "."
                a2 = f[5] if len(f) > 5 and f[5] else a1
                rows.append({"gene": f"HLA-{gene}", "allele_1": a1, "allele_2": a2})
    except Exception as e:
        rows.append({"gene": "error", "allele_1": str(e), "allele_2": ""})

    df = pd.DataFrame(rows)
    with open(a.out, "w") as fh:
        fh.write(f"# hla type for {a.sample} (4-digit, IPD-IMGT/HLA)\n")
    df.to_csv(a.out, sep="\t", index=False, mode="a")
    print(f"[hla_report] {len(df)} genes typed", file=sys.stderr)


if __name__ == "__main__":
    main()
