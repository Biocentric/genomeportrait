#!/usr/bin/env python3
"""Predict ABO, RhD, Kell, Duffy and Kidd blood groups from genotypes.

    python3 blood_groups.py --vcf genotypes.vcf.gz [--rhd rhd_depth.tsv] [--sample S] --out S.blood_groups.tsv
    python3 blood_groups.py --vcf panel.vcf.gz --psam panel.psam --summary-out by_population.tsv   (validation)

The VCF holds the defining variants (assets/blood_groups/sites.pvar), with hom-ref sites
present (from `pgs_engine.py target` on a gVCF). GRCh38, plus strand:

  ABO   9:133257521 T>TC  rs8176719  c.261delG; the reference genome carries the O allele,
                                     so the inserted C marks a functional (A or B) allele
        9:133255935 G>T   rs8176746  c.796C>A, B-specific
        9:133255928 C>G   rs8176747  c.803G>C, B-specific (checks rs8176746)
  KEL   7:142957921 G>A   rs8176058  c.578C>T, K (KEL*01)
  ACKR1 1:159205564 G>A   rs12075    c.125G>A, FY*B (G = FY*A)
        1:159204893 T>C   rs2814778  c.-67T>C, GATA-box: no Duffy on red cells (FY*B background)
  JK    18:45739554 G>A   rs1058396  c.838G>A, JK*B (G = JK*A)

RhD comes from read depth over the RHD gene (--rhd): RhD-negative is, in most people of
European ancestry and many others, a deletion of the whole gene. Two copies read about 1,
one about 0.5, none about 0. RHD variants that silence an intact gene (RHD*Psi, common in
Africa; weak and partial D) are not detected, so 'RhD positive' is the safer statement.

This predicts the common antigens only. It is not a blood-bank typing.
"""
import argparse
import collections
import gzip
import sys

SITES = {
    "abo_o":   ("9", 133257521, "T", "TC"),
    "abo_b1":  ("9", 133255935, "G", "T"),
    "abo_b2":  ("9", 133255928, "C", "G"),
    "kel":     ("7", 142957921, "G", "A"),
    "fy_ab":   ("1", 159205564, "G", "A"),
    "fy_gata": ("1", 159204893, "T", "C"),
    "jk":      ("18", 45739554, "G", "A"),
}


def read_vcf(path):
    """{sample: {site: ALT allele count or None}}"""
    key = {(c, p, r, a): k for k, (c, p, r, a) in SITES.items()}
    out, samples = {}, []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for l in fh:
            if l.startswith("##"):
                continue
            f = l.rstrip("\n").split("\t")
            if l.startswith("#"):
                samples = f[9:]
                out = {s: {} for s in samples}
                continue
            c = f[0][3:] if f[0].startswith("chr") else f[0]
            k = key.get((c, int(f[1]), f[3], f[4]))
            if k is None:
                continue
            for s, g in zip(samples, f[9:]):
                gt = g.split(":")[0].replace("|", "/")
                out[s][k] = None if "." in gt else sum(x == "1" for x in gt.split("/"))
    return out


def abo(g):
    o_fun, b1, b2 = g.get("abo_o"), g.get("abo_b1"), g.get("abo_b2")
    if o_fun is None or b1 is None:
        return "not determined", ""
    n_o = 2 - o_fun                      # O alleles: copies WITHOUT the inserted C
    n_func = o_fun
    note = ""
    if b2 is not None and b2 != b1:
        note = "the two B-defining variants disagree; rs8176746 used"
    n_b = b1
    if n_b > n_func:
        return "not determined", "B variant on an O allele (rare allele)"
    if n_func == 0:
        return "O", note
    if n_b == 0:
        return "A", note
    if n_b == n_func:
        return "B", note
    return "AB", note


def kell(g):
    k = g.get("kel")
    if k is None:
        return "not determined"
    return {0: "K- k+", 1: "K+ k+", 2: "K+ k-"}[k]


def duffy(g):
    ab, gata = g.get("fy_ab"), g.get("fy_gata")
    if ab is None or gata is None:
        return "not determined"
    # GATA -67C silences FY*B on red cells; assume it sits on FY*B alleles (as it does in
    # nearly all carriers) and remove that many FY*B copies before reading the antigens
    fyb = max(0, ab - gata)
    fya = 2 - ab
    if fya == 0 and fyb == 0:
        return "Fy(a- b-)"
    return f"Fy(a{'+' if fya else '-'} b{'+' if fyb else '-'})"


def kidd(g):
    j = g.get("jk")
    if j is None:
        return "not determined"
    return {0: "Jk(a+ b-)", 1: "Jk(a+ b+)", 2: "Jk(a- b+)"}[j]


def rhd_call(path):
    """(phenotype, copies, ratio) from the depth table written by the RHD depth step."""
    if not path:
        return "not determined", "", ""
    vals = {}
    with open(path) as fh:
        for l in fh:
            k, v = l.rstrip("\n").split("\t")[:2]
            vals[k] = v
    ratio = float(vals["ratio"])
    copies = 0 if ratio < 0.2 else 1 if ratio < 0.75 else 2
    return ("RhD negative" if copies == 0 else "RhD positive"), copies, f"{ratio:.2f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--rhd")
    ap.add_argument("--sample")
    ap.add_argument("--out")
    ap.add_argument("--psam", help="validation: population labels for a multi-sample VCF")
    ap.add_argument("--pop-col", default="SuperPop")
    ap.add_argument("--summary-out")
    a = ap.parse_args()
    geno = read_vcf(a.vcf)

    if a.summary_out:
        pop = {}
        with open(a.psam) as fh:
            h = fh.readline().rstrip("\n").lstrip("#").split("\t")
            for l in fh:
                r = l.rstrip("\n").split("\t")
                pop[r[h.index("IID")]] = r[h.index(a.pop_col)]
        tally = collections.defaultdict(collections.Counter)
        for s, g in geno.items():
            p = pop.get(s, "?")
            tally[(p, "ABO")][abo(g)[0]] += 1
            tally[(p, "Kell")][kell(g)] += 1
            tally[(p, "Duffy")][duffy(g)] += 1
            tally[(p, "Kidd")][kidd(g)] += 1
        with open(a.summary_out, "w") as fh:
            fh.write("population\tsystem\tphenotype\tn\tpercent\n")
            for (p, sysname), c in sorted(tally.items()):
                n = sum(c.values())
                for ph, k in sorted(c.items(), key=lambda x: -x[1]):
                    fh.write(f"{p}\t{sysname}\t{ph}\t{k}\t{100 * k / n:.1f}\n")
        return

    if len(geno) != 1:
        sys.exit("blood_groups: expected a single-sample VCF")
    g = next(iter(geno.values()))
    abo_ph, abo_note = abo(g)
    rh, copies, ratio = rhd_call(a.rhd)
    rows = [
        ("ABO", abo_ph, abo_note or "from rs8176719 (O), rs8176746 and rs8176747 (B)"),
        ("RhD", rh, f"RHD gene copies from read depth: {copies} (depth ratio {ratio}); "
                    "silencing variants in an intact RHD gene are not detected" if ratio else
                    "no read-depth estimate"),
        ("Kell", kell(g), "rs8176058 (K/k)"),
        ("Duffy", duffy(g), "rs12075 (Fya/Fyb), rs2814778 (GATA-box, Fy(b) silenced on red cells)"),
        ("Kidd", kidd(g), "rs1058396 (Jka/Jkb)"),
    ]
    with open(a.out, "w") as fh:
        fh.write(f"# blood groups for {a.sample} - predicted from genotypes, not a blood-bank typing\n")
        fh.write("system\tpredicted\tbasis\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")


if __name__ == "__main__":
    main()
