#!/usr/bin/env python3
"""
trait_lookup.py — genotype the curated non-clinical trait panel in a sample VCF.

Pure-Python: reads the (optionally bgzipped) VCF directly — no bcftools needed — so it
runs in the slim report container. For each panel rsID we find the record by its ID
(or by Existing_variation from VEP), read the sample genotype, count effect-allele
dosage, and emit the matching plain-language interpretation from the panel.
"""
import argparse
import gzip
import sys

import pandas as pd


def opener(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") or str(p).endswith(".bgz") else open(p)


def load_genotypes(vcf, rsids):
    """Return {rsid: set(called_alleles)} for panel rsIDs present in the VCF."""
    want = set(rsids)
    found = {}
    with opener(vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 10:
                continue
            chrom, pos, vid, ref, alt = f[0], f[1], f[2], f[3], f[4]
            ids = set(vid.split(";"))
            # also scan VEP CSQ Existing_variation in INFO for rsIDs
            if "rs" in f[7]:
                for tok in f[7].replace("|", ";").replace(",", ";").split(";"):
                    if tok.startswith("rs"):
                        ids.add(tok)
            hit = ids & want
            if not hit:
                continue
            alleles = [ref] + alt.split(",")
            fmt = f[8].split(":")
            try:
                gt_idx = fmt.index("GT")
            except ValueError:
                continue
            gt = f[9].split(":")[gt_idx].replace("|", "/")
            called = set()
            for tok in gt.split("/"):
                if tok in (".", ""):
                    continue
                try:
                    called.add(alleles[int(tok)])
                except (ValueError, IndexError):
                    pass
            for rs in hit:
                found[rs] = called
    return found


def classify(eff, oth, called):
    if not called:
        return "no_call", "./."
    if called == {eff}:
        return "hom_effect", f"{eff}/{eff}"
    if called == {oth}:
        return "hom_other", f"{oth}/{oth}"
    if eff in called and oth in called:
        return "het", f"{eff}/{oth}"
    return "other_allele", "/".join(sorted(called))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    panel = pd.read_csv(a.panel, sep="\t", comment="#")
    gts = load_genotypes(a.vcf, panel["rsid"].tolist())

    rows = []
    for _, r in panel.iterrows():
        rs = r["rsid"]
        called = gts.get(rs, set())
        state, dosage = classify(str(r["effect_allele"]), str(r["other_allele"]), called)
        interp = {
            "hom_effect": r["geno_hom_effect"],
            "het": r["geno_het"],
            "hom_other": r["geno_hom_other"],
            "other_allele": "Genotype outside modelled alleles — see source",
            "no_call": "Not covered / no confident call",
        }[state]
        rows.append({
            "category": r["category"], "trait": r["trait"], "gene": r["gene"],
            "rsid": rs, "your_genotype": dosage, "interpretation": interp,
            "effect_allele": r["effect_allele"],
        })

    df = pd.DataFrame(rows).sort_values(["category", "trait"])
    with open(a.out, "w") as fh:
        fh.write(f"# trait_profile for {a.sample}\n")
    df.to_csv(a.out, sep="\t", index=False, mode="a")
    called_n = (df["your_genotype"] != "./.").sum()
    print(f"[trait_lookup] {called_n}/{len(df)} panel SNPs genotyped", file=sys.stderr)


if __name__ == "__main__":
    main()
