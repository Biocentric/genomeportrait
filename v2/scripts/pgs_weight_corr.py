#!/usr/bin/env python3
"""Usage: python3 pgs_weight_corr.py <scorefile_dir> ANCHOR:TEST1,TEST2 [ANCHOR:...]

Correlate effect weights between scoring files on shared, allele-aligned variants.
A negative correlation against an anchor of known direction = reversed trait coding."""
import gzip, sys, math, os
D = sys.argv[1]
COMP = str.maketrans("ACGT", "TGCA")
def load(p):
    out = {}
    with gzip.open(os.path.join(D, p + ".txt.gz"), "rt") as fh:
        hdr = None
        for line in fh:
            if line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")
            if hdr is None: hdr = {k: i for i, k in enumerate(f)}; continue
            try:
                c, pos = f[hdr["hm_chr"]], f[hdr["hm_pos"]]
                ea = f[hdr["effect_allele"]].upper()
                oa = (f[hdr["other_allele"]] if "other_allele" in hdr else "") or \
                     (f[hdr["hm_inferOtherAllele"]] if "hm_inferOtherAllele" in hdr else "")
                w = float(f[hdr["effect_weight"]])
            except (KeyError, ValueError, IndexError):
                continue
            if not c or not pos: continue
            out[(c, pos)] = (ea, oa.upper(), w)
    return out
def corr(a, b):
    xs, ys = [], []
    for k, (ea1, oa1, w1) in a.items():
        if k not in b: continue
        ea2, oa2, w2 = b[k]
        if ea1 == ea2 or ea1 == ea2.translate(COMP): s = 1
        elif ea1 == oa2 or ea1 == oa2.translate(COMP) or (oa1 and oa1 == ea2): s = -1
        else: continue
        xs.append(w1); ys.append(s * w2)
    n = len(xs)
    if n < 30: return n, float("nan")
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs)); sy = math.sqrt(sum((y - my) ** 2 for y in ys))
    return n, sxy / (sx * sy)
cache = {}
def L(p):
    if p not in cache: cache[p] = load(p)
    return cache[p]
for line in sys.argv[2:]:
    anchor, tests = line.split(":")
    print(f"anchor {anchor}")
    for t in tests.split(","):
        n, r = corr(L(anchor), L(t))
        print(f"   {t}: shared={n:>7}  r={r:+.3f}")
