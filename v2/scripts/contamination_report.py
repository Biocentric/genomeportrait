#!/usr/bin/env python3
"""VerifyBamID2 result as genome-facts rows.

    python3 contamination_report.py --sample S --selfsm S.selfSM --out S.contamination_facts.tsv

FREEMIX is the estimated fraction of reads from another person. Below 2% does not affect
genotypes at normal depth; 2-5% can create false heterozygous calls at low allele
fraction; above 5% the genotypes and everything built on them are unreliable.
"""
import argparse
import csv
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--selfsm", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    with open(a.selfsm) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        sys.exit("contamination_report: empty selfSM")
    r = {k.lstrip("#"): v for k, v in rows[0].items()}
    freemix = float(r["FREEMIX"])
    if freemix < 0.02:
        verdict = "clean: no effect on genotypes expected"
    elif freemix < 0.05:
        verdict = "check: can add false heterozygous calls at low allele fraction"
    else:
        verdict = "contaminated: genotypes and everything derived from them are unreliable"
    with open(a.out, "w") as fh:
        fh.write(f"# contamination for {a.sample}\n")
        fh.write("fact\tvalue\tnote\n")
        fh.write(f"DNA from another person (FREEMIX)\t{100 * freemix:.2f}%\t"
                 f"{verdict}; VerifyBamID2 over {int(float(r['SNPS'])):,} common SNPs at {float(r['AVG_DP']):.1f}x\n")
    print(f"[contamination] {a.sample}: FREEMIX {freemix:.4f}", file=sys.stderr)


if __name__ == "__main__":
    main()
