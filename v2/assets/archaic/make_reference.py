#!/usr/bin/env python3
"""Per-population archaic sequence from the published hmmix calls for 1000 Genomes.

    python3 assets/archaic/make_reference.py hg38_1000g_segments.txt > assets/archaic/1000g_reference.tsv

Input: Skov's hmmix segments for the 1000 Genomes high-coverage genomes (Zenodo 14136628,
CC-BY 4.0), one row per segment and phased haplotype. Per person, the segments of both
haplotypes with mean posterior >= 0.8 are merged, which is what an unphased (diploid) run
reports: sequence where at least one copy is archaic. The table holds the median, 10th and
90th percentile per population. African populations are absent: they are hmmix's outgroup.
"""
import collections
import statistics
import sys

MIN_PROB = 0.8
segs = collections.defaultdict(list)
pop = {}
with open(sys.argv[1]) as fh:
    h = fh.readline().rstrip("\n").split("\t")
    ix = {k: h.index(k) for k in ("name", "pop", "region", "chrom", "start", "end", "mean_prob")}
    for l in fh:
        f = l.rstrip("\n").split("\t")
        if float(f[ix["mean_prob"]]) < MIN_PROB:
            continue
        segs[f[ix["name"]]].append((f[ix["chrom"]], int(f[ix["start"]]), int(f[ix["end"]])))
        pop[f[ix["name"]]] = (f[ix["pop"]], f[ix["region"]])


def union_bp(iv):
    total, cur = 0, None
    for c, a, b in sorted(iv):
        if cur and cur[0] == c and a <= cur[2]:
            cur[2] = max(cur[2], b)
        else:
            if cur:
                total += cur[2] - cur[1]
            cur = [c, a, b]
    return total + (cur[2] - cur[1] if cur else 0)


by_pop = collections.defaultdict(list)
for s, iv in segs.items():
    by_pop[pop[s]].append(union_bp(iv) / 1e6)

print("# archaic sequence (Mb, at least one copy, mean posterior >= 0.8) per 1000 Genomes population;")
print("# from Skov, hmmix calls for 1000 Genomes hg38, Zenodo 14136628 (CC-BY 4.0)")
print("population\tregion\tn\tmedian_Mb\tp10_Mb\tp90_Mb")
for (p, r), v in sorted(by_pop.items(), key=lambda x: (x[0][1], x[0][0])):
    q = statistics.quantiles(v, n=10)
    print(f"{p}\t{r}\t{len(v)}\t{statistics.median(v):.1f}\t{q[0]:.1f}\t{q[-1]:.1f}")
