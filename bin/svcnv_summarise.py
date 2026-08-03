#!/usr/bin/env python3
"""svcnv_summarise.py — compact counts of structural and copy-number variants."""
import argparse
import gzip
import sys

import pandas as pd


def opener(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def sv_counts(path):
    counts = {}
    if not path:
        return counts
    try:
        with opener(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                info = line.split("\t")[7] if len(line.split("\t")) > 7 else ""
                kv = dict(x.split("=", 1) for x in info.split(";") if "=" in x)
                svt = kv.get("SVTYPE", "OTHER")
                counts[svt] = counts.get(svt, 0) + 1
    except Exception as e:
        print(f"[svcnv] SV parse: {e}", file=sys.stderr)
    return counts


def cnv_counts(path):
    """Parse a CNVpytor call TSV (no header): col0 = type (deletion|duplication)."""
    if not path:
        return (0, 0)
    gains = losses = 0
    try:
        with open(path) as fh:
            for line in fh:
                t = line.split("\t", 1)[0].strip().lower()
                if t.startswith("dup"):
                    gains += 1
                elif t.startswith("del"):
                    losses += 1
    except Exception as e:
        print(f"[svcnv] CNV parse: {e}", file=sys.stderr)
    return (gains, losses)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--sv", default=None)
    ap.add_argument("--cnr", default=None)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    rows = []
    for svt, n in sorted(sv_counts(a.sv).items()):
        rows.append({"category": "Structural variant", "type": svt, "count": n})
    dup, dele = cnv_counts(a.cnr)
    rows.append({"category": "Copy-number (CNVpytor)", "type": "duplications", "count": dup})
    rows.append({"category": "Copy-number (CNVpytor)", "type": "deletions", "count": dele})
    if not rows:
        rows.append({"category": "none", "type": "none", "count": 0})

    df = pd.DataFrame(rows)
    with open(a.out, "w") as fh:
        fh.write(f"# svcnv_summary for {a.sample}\n")
    df.to_csv(a.out, sep="\t", index=False, mode="a")
    print("[svcnv_summarise] done", file=sys.stderr)


if __name__ == "__main__":
    main()
