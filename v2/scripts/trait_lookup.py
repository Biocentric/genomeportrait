#!/usr/bin/env python3
"""
trait_lookup.py — interpret the curated non-clinical trait panel from targeted genotypes.

Reads a per-position genotype table (CHROM POS REF ALT GT DP) produced by bcftools mpileup
at the panel coordinates. Because the sites are genotyped directly from the BAM (not taken
from a variant-only VCF), homozygous-reference calls are captured too — so a panel SNP is
reported as hom-reference rather than "no call" when the individual simply matches the genome.
"""
import argparse
import sys

import pandas as pd


def load_geno(path):
    """Return {(chrom,pos): set(called_alleles)} from the bcftools query TSV."""
    g = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            # f layout: CHROM POS REF ALT GT [DP]
            chrom, pos, ref, alt, gt = f[0], f[1], f[2], f[3], f[4]
            alleles = [ref] + ([] if alt in (".", "") else alt.split(","))
            called = set()
            for tok in gt.replace("|", "/").split("/"):
                if tok in (".", ""):
                    continue
                try:
                    called.add(alleles[int(tok)])
                except (ValueError, IndexError):
                    pass
            g[(chrom, str(pos))] = called
    return g


def _comp(b):
    return {"A": "T", "T": "A", "C": "G", "G": "C"}.get(b, b)


def classify(eff, oth, called):
    if not called:
        return "no_call", "./."
    dosage = "/".join(sorted(called))
    # Match on the forward strand, or the complement (panel alleles are sometimes given
    # on the opposite strand to the GATK hg38 forward-strand VCF, e.g. rs4988235 T/C vs A/G).
    for e, o in [(eff, oth), (_comp(eff), _comp(oth))]:
        if called == {e}:
            return "hom_effect", dosage
        if called == {o}:
            return "hom_other", dosage
        if e in called and o in called:
            return "het", dosage
    return "other_allele", dosage


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--geno", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    panel = pd.read_csv(a.panel, sep="\t", comment="#")
    geno = load_geno(a.geno)

    rows = []
    for _, r in panel.iterrows():
        key = (str(r["chrom"]), str(r["pos_grch38"]))
        called = geno.get(key, set())
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
            "rsid": r["rsid"], "your_genotype": dosage, "interpretation": interp,
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
