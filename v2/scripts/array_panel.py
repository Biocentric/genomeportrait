#!/usr/bin/env python3
"""Site table for a 23andMe-format export, from the Illumina GSA-24 v3.0 manifests.

    python3 array_panel.py --b37 GSA-24v3-0_A1.csv --b38 GSA-24v3-0_A2.csv --out gsa_v3_sites.tsv.gz

The GSA is the array 23andMe's v5 chip is built on, so its SNPs are what downstream tools
that read 23andMe files expect. The A1 manifest gives GRCh37 positions (the 23andMe
format's build), A2 the GRCh38 positions this pipeline genotypes at. A probe whose
RefStrand differs between the two builds sits in a region that changed orientation; its
genotype is complemented on export. Indel probes, probes on different chromosomes in the
two builds, and probes without a position are left out.

Output columns: name chrom pos37 pos38 flip (chrom: 1-22, X, Y, MT; XY probes are X).
"""
import argparse
import csv
import gzip
import sys

CHROMS = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]


def load(path):
    out = {}
    with open(path, newline="") as fh:
        for line in fh:
            if line.startswith("[Assay]"):
                break
        r = csv.reader(fh)
        h = next(r)
        i = {k: n for n, k in enumerate(h)}
        for row in r:
            if not row or row[0].startswith("[Controls]"):
                break
            out[row[i["IlmnID"]]] = (row[i["Name"]], row[i["SNP"]], row[i["Chr"]], row[i["MapInfo"]],
                                     row[i["RefStrand"]])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--b37", required=True)
    ap.add_argument("--b38", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    b37, b38 = load(a.b37), load(a.b38)
    rows, skipped = {}, dict(indel=0, chrom=0, position=0, duplicate=0)
    for k in sorted(b37.keys() & b38.keys()):      # sorted: the same table on every run
        name, snp, c37, p37, s37 = b37[k]
        _n, _s, c38, p38, s38 = b38[k]
        if snp in ("[I/D]", "[D/I]"):
            skipped["indel"] += 1
            continue
        c37 = "X" if c37 == "XY" else c37
        c38 = "X" if c38 == "XY" else c38
        if c37 != c38 or c37 not in CHROMS:
            skipped["chrom"] += 1
            continue
        if not p37.isdigit() or not p38.isdigit() or p37 == "0" or p38 == "0":
            skipped["position"] += 1
            continue
        key = (c38, int(p38))
        if key in rows:
            # several probes for one SNP: keep one, preferring an rs name (what tools that
            # read 23andMe files look up), then the alphabetically first
            skipped["duplicate"] += 1
            old = rows[key][0]
            if (not name.startswith("rs"), name) >= (not old.startswith("rs"), old):
                continue
        flip = int(s37 in "+-" and s38 in "+-" and s37 != s38)
        rows[key] = (name, c37, int(p37), int(p38), flip)
    order = {c: n for n, c in enumerate(CHROMS)}
    with gzip.open(a.out, "wt") as fh:
        fh.write("name\tchrom\tpos37\tpos38\tflip\n")
        for name, c, p37, p38, flip in sorted(rows.values(), key=lambda r: (order[r[1]], r[3])):
            fh.write(f"{name}\t{c}\t{p37}\t{p38}\t{flip}\n")
    print(f"[array_panel] {len(rows):,} SNP sites; skipped " +
          ", ".join(f"{k} {v:,}" for k, v in skipped.items()), file=sys.stderr)


if __name__ == "__main__":
    main()
