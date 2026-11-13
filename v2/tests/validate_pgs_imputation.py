#!/usr/bin/env python3
"""How well do pooled and population-specific frequencies stand in for missing genotypes?

    python3 validate_pgs_imputation.py --ref-prefix ref/ref --scores scores --target target.vcf.gz \
        --pooled ref/ref.afreq --pop-afreq-dir popfreq --ref-scores ref_scores \
        --plink2 "singularity exec plink2.sif plink2" --workdir work

The sites a real sample is missing (read from its target VCF) are hidden in every
reference genome. For each reference genome and score, the true contribution of those
sites is compared with what filling them in would give: 2 x the pooled ALT frequency, or
2 x the frequency of the genome's own super-population. Errors are expressed in SDs of the
score within that population, so they read like a shift in Z.

The genome's own genotypes are part of its population's frequencies (~700 genomes per
population), which flatters the population-specific estimate very slightly.
"""
import argparse
import glob
import gzip
import os
import subprocess
import sys

import numpy as np


def read_afreq(path):
    with open(path) as fh:
        h = fh.readline().rstrip("\n").lstrip("#").split("\t")
        i = h.index("ALT_FREQS")
        return np.array([float(l.split("\t")[i]) if l.split("\t")[i] not in ("NA", "nan") else 0.0 for l in fh])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref-prefix", required=True)
    ap.add_argument("--scores", required=True)
    ap.add_argument("--target", required=True)
    ap.add_argument("--pooled", required=True)
    ap.add_argument("--pop-afreq-dir", required=True)
    ap.add_argument("--ref-scores", required=True)
    ap.add_argument("--plink2", required=True)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--pop-col", default="SuperPop")
    a = ap.parse_args()
    os.makedirs(a.workdir, exist_ok=True)

    ids, alts = [], []
    with open(a.ref_prefix + ".pvar") as fh:
        for l in fh:
            if l[0] != "#":
                f = l.split("\t", 5)
                ids.append(f[2]); alts.append(f[4].rstrip("\n"))
    n = len(ids)
    missing = np.zeros(n, dtype=bool)
    with gzip.open(a.target, "rt") as fh:
        k = 0
        for l in fh:
            if l[0] == "#":
                continue
            if l.rstrip("\n").endswith("./."):
                missing[k] = True
            k += 1
    if k != n:
        sys.exit(f"target has {k} variants, reference {n}")
    rows = np.flatnonzero(missing)
    pos_of = {r: j for j, r in enumerate(rows)}
    print(f"{len(rows):,} of {n:,} sites missing in the target are hidden in every reference genome", file=sys.stderr)

    # one weight matrix over the hidden sites, one column per score
    files = sorted(glob.glob(os.path.join(a.scores, "*.npz")))
    pids, W = [], np.zeros((len(rows), len(files)))
    for j, p in enumerate(files):
        z = np.load(p)
        m = missing[z["row"]]
        for r, w in zip(z["row"][m], z["w"][m]):
            W[pos_of[r], j] += w
        pids.append(os.path.basename(p)[:-4])
    keep = np.flatnonzero(np.abs(W).sum(0) > 0)
    pids = [pids[j] for j in keep]
    W = W[:, keep]
    sf = os.path.join(a.workdir, "hidden.score")
    with open(sf, "w") as fh:
        fh.write("ID\tA1\t" + "\t".join(pids) + "\n")
        for j, r in enumerate(rows):
            fh.write(f"{ids[r]}\t{alts[r]}\t" + "\t".join(f"{x:.10g}" for x in W[j]) + "\n")
    out = os.path.join(a.workdir, "hidden")
    cmd = (f"{a.plink2} --threads 16 --pfile {a.ref_prefix} --read-freq {a.pooled} "
           f"--score {sf} 1 2 header-read cols=scoresums --score-col-nums 3-{len(pids) + 2} --out {out}")
    subprocess.run(cmd, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    with open(out + ".sscore") as fh:
        h = fh.readline().rstrip("\n").lstrip("#").split("\t")
        cols = [h.index(f"{p}_SUM") for p in pids]
        iid = h.index("IID")
        true = {}
        for l in fh:
            f = l.rstrip("\n").split("\t")
            true[f[iid]] = np.array([float(f[c]) for c in cols])

    pop = {}
    with open(a.ref_prefix + ".psam") as fh:
        h = fh.readline().rstrip("\n").lstrip("#").split("\t")
        for l in fh:
            f = l.rstrip("\n").split("\t")
            pop[f[h.index("IID")]] = f[h.index(a.pop_col)]
    pooled = read_afreq(a.pooled)[rows]
    exp_pooled = 2 * pooled @ W
    pops = sorted(set(pop.values()))
    print(f"\n{'population':10s} {'n':>5s}  {'bias pooled':>11s} {'bias own':>9s}  {'RMSE pooled':>11s} {'RMSE own':>9s}  "
          f"{'own better':>10s}   (median over {len(pids)} scores, in SD of the score within the population)")
    for p in pops:
        members = [i for i in true if pop.get(i) == p]
        T = np.array([true[i] for i in members])
        exp_own = 2 * read_afreq(os.path.join(a.pop_afreq_dir, f"{p}.afreq"))[rows] @ W
        # SD of each full score within the population, from the reference scoring
        sd = []
        for pid in pids:
            with open(os.path.join(a.ref_scores, f"{pid}.sscore")) as fh:
                h = fh.readline().rstrip("\n").lstrip("#").split("\t")
                ci, cs = h.index("IID"), next(i for i, c in enumerate(h) if c.endswith("_SUM") and "DOSAGE" not in c)
                v = [float(f[cs]) for f in (l.rstrip("\n").split("\t") for l in fh) if pop.get(f[ci]) == p]
            sd.append(np.std(v))
        sd = np.array(sd)
        ok = sd > 0
        e_pool = (T - exp_pooled)[:, ok] / sd[ok]
        e_own = (T - exp_own)[:, ok] / sd[ok]
        bias_pool, bias_own = np.abs(e_pool.mean(0)), np.abs(e_own.mean(0))
        rmse_pool, rmse_own = np.sqrt((e_pool ** 2).mean(0)), np.sqrt((e_own ** 2).mean(0))
        better = np.mean(rmse_own <= rmse_pool)
        print(f"{p:10s} {len(members):5d}  {np.median(bias_pool):11.4f} {np.median(bias_own):9.4f}  "
              f"{np.median(rmse_pool):11.4f} {np.median(rmse_own):9.4f}  {100 * better:9.0f}%")
        print(f"{'':10s} {'':5s}  max bias pooled {bias_pool.max():.3f}, own {bias_own.max():.3f}; "
              f"p95 RMSE pooled {np.percentile(rmse_pool, 95):.3f}, own {np.percentile(rmse_own, 95):.3f}")


if __name__ == "__main__":
    main()
