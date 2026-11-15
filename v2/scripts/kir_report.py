#!/usr/bin/env python3
"""KIR gene content and alleles from T1K, with the A/B haplotype classification.

    python3 kir_report.py --sample S --t1k S_genotype.tsv --out S.kir.tsv

KIR genes code for receptors on natural killer cells. People differ in which genes they
carry, not only in alleles. Haplotypes come in two groups: A haplotypes carry a fixed set of
mostly inhibitory genes; B haplotypes add one or more of KIR2DL2, KIR2DL5, KIR2DS1, KIR2DS2,
KIR2DS3, KIR2DS5 and KIR3DS1. A genome with none of those is A/A; otherwise it carries at
least one B haplotype (B/x). The framework genes (KIR3DL3, KIR3DP1, KIR2DL4, KIR3DL2) are
present in almost everyone; if any is missing, the typing is flagged.

A gene with one allele listed carries one copy or two identical ones. KIR2DL2/KIR2DL3 and
KIR3DL1/KIR3DS1 are each alleles of one locus, so a genome with both of a pair has one of each.
"""
import argparse
import sys

FRAMEWORK = ["KIR3DL3", "KIR3DP1", "KIR2DL4", "KIR3DL2"]
B_GENES = ["KIR2DL2", "KIR2DL5A", "KIR2DL5B", "KIR2DS1", "KIR2DS2", "KIR2DS3", "KIR2DS5", "KIR3DS1"]
ORDER = ["KIR3DL3", "KIR2DS2", "KIR2DL2", "KIR2DL3", "KIR2DL5B", "KIR2DS3", "KIR2DP1", "KIR2DL1",
         "KIR3DP1", "KIR2DL4", "KIR3DL1", "KIR3DS1", "KIR2DL5A", "KIR2DS5", "KIR2DS1", "KIR2DS4", "KIR3DL2"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--t1k", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    calls = {}
    with open(a.t1k) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            # T1K: gene, number of alleles, allele1, abundance1, quality1, allele2, ...
            # Quality 0 marks a call made from a handful of paralogous reads: not counted.
            alleles = [f[i] for i in (2, 5)
                       if len(f) > i + 2 and f[i] not in (".", "") and int(f[i + 2]) > 0]
            if int(f[1]) > 0 and alleles:
                calls[f[0]] = alleles

    present = [g for g in ORDER if g in calls]
    missing_frame = [g for g in FRAMEWORK if g not in calls]
    b_present = [g for g in B_GENES if g in calls]
    haplo = "A/A" if not b_present else "B/x (at least one B haplotype)"

    with open(a.out, "w") as fh:
        fh.write(f"# KIR genes for {a.sample} (T1K, IPD-KIR); not a clinical typing\n")
        fh.write("gene\tpresent\tallele_1\tallele_2\n")
        for g in ORDER:
            al = calls.get(g, [])
            # one allele reported: one copy, or two identical ones (depth cannot tell them
            # apart reliably), so the second column stays empty rather than repeating it
            fh.write(f"{g}\t{'yes' if al else 'no'}\t{al[0] if al else ''}\t"
                     f"{al[1] if len(al) > 1 else ''}\n")
        fh.write(f"haplotype group\t{haplo}\t{', '.join(b_present) or 'no B-haplotype genes'}\t\n")
        if missing_frame:
            fh.write(f"check\tframework gene(s) not found: {', '.join(missing_frame)}\t"
                     "low coverage or reads outside the extracted region\t\n")
    print(f"[kir_report] {a.sample}: {len(present)} genes, {haplo}"
          + (f"; framework missing {missing_frame}" if missing_frame else ""), file=sys.stderr)


if __name__ == "__main__":
    main()
