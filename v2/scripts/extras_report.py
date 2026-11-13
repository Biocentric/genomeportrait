#!/usr/bin/env python3
"""extras_report.py — combine telomere-length estimate and mtDNA heteroplasmy count."""
import argparse
import sys

import pandas as pd


def telomere_len(path):
    if not path:
        return None
    try:
        df = pd.read_csv(path, sep="\t", comment="#")
        col = next((c for c in df.columns if "LENGTH" in c.upper() or "tel" in c.lower()), None)
        if col:
            return round(float(df[col].dropna().iloc[0]), 2)
    except Exception as e:
        print(f"[extras] telseq parse: {e}", file=sys.stderr)
    return None


def heteroplasmy_count(path):
    if not path:
        return None
    try:
        df = pd.read_csv(path, sep="\t")
        lvl = next((c for c in df.columns if "level" in c.lower()), None)
        if lvl:
            sub = df[(df[lvl].astype(float) >= 0.01) & (df[lvl].astype(float) <= 0.95)]
            return int(len(sub))
        return int(len(df))
    except Exception as e:
        print(f"[extras] mutserve parse: {e}", file=sys.stderr)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--telseq", default=None)
    ap.add_argument("--mito", default=None)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    rows = []
    tl = telomere_len(a.telseq)
    if tl is not None:
        rows.append({
            "metric": "Mean telomere length",
            "value": f"{tl} kb",
            "what_it_measures":
                "Average length of the repetitive TTAGGG caps across all chromosome ends, "
                "estimated from how many telomeric-repeat reads are in the whole-genome data "
                "(Telseq). It is one genome-wide average, not a per-chromosome measurement.",
            "how_to_read_it":
                "Telseq values for adult blood WGS typically land roughly in the 3-6 kb range, "
                "and shorter averages are seen at older ages. The absolute number is strongly "
                "method-dependent: Telseq, qPCR and Southern blot give different values for the "
                "same sample, so it is only comparable with other Telseq WGS estimates produced "
                "the same way — not against a number from another lab or method.",
            "tool": "Telseq",
        })
    hc = heteroplasmy_count(a.mito)
    if hc is not None:
        rows.append({
            "metric": "mtDNA heteroplasmies (1-95%)",
            "value": hc,
            "what_it_measures":
                "Positions where the mitochondrial genome is NOT uniform: cells carry a mix of "
                "two mtDNA sequences. Each cell holds hundreds of mitochondrial genomes, so a "
                "variant can be present in only a fraction of them.",
            "how_to_read_it":
                "0 means no position had a minor allele between 1% and 95% at sufficient depth — "
                "i.e. the mitochondrial genome looks homoplasmic (uniform). Low counts are the "
                "common finding; detection depends on mtDNA coverage.",
            "tool": "mutserve",
        })
    if not rows:
        rows.append({"metric": "extras", "value": "none",
                     "what_it_measures": "modules skipped or no data",
                     "how_to_read_it": ".", "tool": "."})

    df = pd.DataFrame(rows)
    with open(a.out, "w") as fh:
        fh.write(f"# extras for {a.sample}\n")
    df.to_csv(a.out, sep="\t", index=False, mode="a")
    print("[extras_report] done", file=sys.stderr)


if __name__ == "__main__":
    main()
