#!/usr/bin/env python3
"""str_report.py — produce a CODIS-style STR genotype table + repeat-expansion sizes."""
import argparse
import gzip
import sys

import pandas as pd

# The 20 expanded CODIS core loci (markers used in forensic identity databases)
CODIS = {"CSF1PO", "FGA", "TH01", "TPOX", "VWA", "D3S1358", "D5S818", "D7S820",
         "D8S1179", "D13S317", "D16S539", "D18S51", "D21S11", "D1S1656",
         "D2S441", "D2S1338", "D10S1248", "D12S391", "D19S433", "D22S1045"}


def opener(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def parse_vcf_str(path, kind):
    rows = []
    try:
        with opener(path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 10:
                    continue
                info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
                fmt = f[8].split(":")
                smp = dict(zip(fmt, f[9].split(":")))
                locus = info.get("GENE") or f[2]
                if kind == "codis":
                    rows.append({"marker": locus, "chrom": f[0], "pos": f[1],
                                 "GB": smp.get("GB", smp.get("GT", "."))})
                else:
                    rows.append({"locus": locus, "repeat_unit": info.get("RU", "."),
                                 "genotype_repeats": smp.get("REPCN", smp.get("GT", "."))})
    except Exception as e:
        print(f"[str_report] {kind} parse error: {e}", file=sys.stderr)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--hipstr", required=True)
    ap.add_argument("--expansionhunter", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    codis = [r for r in parse_vcf_str(a.hipstr, "codis")
             if any(c in str(r["marker"]).upper() for c in CODIS)] or parse_vcf_str(a.hipstr, "codis")
    eh = parse_vcf_str(a.expansionhunter, "eh")

    with open(a.out, "w") as fh:
        fh.write(f"# str_profile for {a.sample}\n")
        fh.write("## CODIS core STR profile (forensic identity markers)\n")
        pd.DataFrame(codis).to_csv(fh, sep="\t", index=False)
        fh.write("\n## Repeat expansion loci (ExpansionHunter)\n")
        pd.DataFrame(eh).to_csv(fh, sep="\t", index=False)
    print(f"[str_report] {len(codis)} STR + {len(eh)} expansion loci", file=sys.stderr)


if __name__ == "__main__":
    main()
