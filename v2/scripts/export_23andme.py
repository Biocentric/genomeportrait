#!/usr/bin/env python3
"""Write a genome as a 23andMe-format raw data file, at the GSA array's SNPs.

    python3 export_23andme.py --gvcf S.g.vcf.gz --sites gsa_v3_sites.tsv.gz --fasta GRCh38.fa \
        --sex XY --sample S --threads 8 --out S.23andme.txt

Genotypes come from the gVCF, so a site where the genome matches the reference is written
as the reference base twice (from the FASTA), not left out; a site the gVCF has no reads
for is "--". Positions are GRCh37 (the format's build), bases are plus strand. For XY
samples chrX outside the pseudo-autosomal regions, chrY and chrM are haploid (one letter);
a heterozygous call there is written "--". Indels are not exported.

The file is for personal use with tools that read 23andMe raw data. It is not a 23andMe
file and says so in its header.
"""
import argparse
import datetime
import gzip
import multiprocessing as mp
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pgs_engine import bgzf_lines, read_tbi  # noqa: E402  (staged next to this script)

COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}
PAR38 = [(10001, 2781479), (155701383, 156030895)]      # chrX pseudo-autosomal, GRCh38
OUT_ORDER = [str(i) for i in range(1, 23)] + ["X", "Y", "MT"]


class Fasta:
    """Single-base lookups through the .fai index."""
    def __init__(self, path):
        self.fh = open(path, "rb")
        self.idx = {}
        with open(path + ".fai") as f:
            for l in f:
                n, ln, off, lb, lw = l.split("\t")[:5]
                self.idx[n] = (int(ln), int(off), int(lb), int(lw))

    def base(self, chrom, pos):
        ln, off, lb, lw = self.idx[chrom]
        p = pos - 1
        self.fh.seek(off + (p // lb) * lw + p % lb)
        return self.fh.read(1).decode().upper()


SITES = {}   # gvcf contig -> list of (pos38, row index)
OPTS = {}


def genotype_chrom(job):
    contig, voffset = job
    sites = SITES[contig]
    fa = Fasta(OPTS["fasta"])
    out = {}
    i, n = 0, len(sites)
    lines = bgzf_lines(OPTS["gvcf"], voffset) if voffset is not None else iter(())
    for raw in lines:
        f = raw.split(b"\t", 9)
        if f[0].decode() != contig:
            break
        pos = int(f[1])
        while i < n and sites[i][0] < pos:
            i += 1                               # not covered by any record: stays missing
        if i >= n:
            break
        ref = f[3].decode().upper()
        info = f[7]
        end = pos + len(ref) - 1
        k = info.find(b"END=")
        if k != -1:
            j = info.find(b";", k)
            end = int(info[k + 4:j if j != -1 else len(info)])
        alts = [x for x in f[4].decode().upper().split(",") if x not in ("<*>", "<NON_REF>", ".")]
        fmt = f[8].split(b":")
        smp = f[9].rstrip(b"\n").split(b":")
        gt = smp[0].decode().replace("|", "/")
        dp = -1
        for key, val in zip(fmt, smp):
            if key in (b"MIN_DP", b"DP"):
                dp = int(val) if val.isdigit() else -1
        if "." in gt:
            continue
        try:
            alleles = [int(x) for x in gt.split("/")]
        except ValueError:
            continue
        if not alts:                             # reference block
            if dp == 0 or any(alleles):
                continue
            while i < n and sites[i][0] <= end:
                b = fa.base(contig, sites[i][0])
                if b in COMP:
                    out[sites[i][1]] = b + b
                i += 1
        else:                                    # variant record: only its own position
            while i < n and sites[i][0] == pos:
                bases = []
                for x in alleles:
                    if x == 0:
                        bases.append(ref[0])
                    elif x <= len(alts):
                        alt = alts[x - 1]
                        if len(alt) == 1 and len(ref) == 1:
                            bases.append(alt)
                        elif alt[:1] == ref[:1]:
                            bases.append(ref[0])     # indel anchored on an unchanged base
                        else:
                            bases = None
                            break
                    else:
                        bases = None
                        break
                if bases and all(b in COMP for b in bases):
                    out[sites[i][1]] = "".join(sorted(bases))
                i += 1
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gvcf", required=True)
    ap.add_argument("--sites", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--sex", default="unknown", help="XX, XY or unknown (unknown: diploid everywhere)")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    starts = read_tbi(a.gvcf + ".tbi")
    chr_style = any(c.startswith("chr") for c in starts)
    contig_of = {c: (("chr" + c) if chr_style else c) for c in OUT_ORDER}
    contig_of["MT"] = "chrM" if chr_style else "MT"
    rows = []
    with gzip.open(a.sites, "rt") as fh:
        next(fh)
        for k, l in enumerate(fh):
            name, c, p37, p38, flip = l.rstrip("\n").split("\t")
            rows.append((name, c, int(p37), int(p38), flip == "1"))
            SITES.setdefault(contig_of[c], []).append((int(p38), k))
    for v in SITES.values():
        v.sort()
    OPTS.update(gvcf=a.gvcf, fasta=a.fasta)
    jobs = [(c, starts.get(c)) for c in SITES]
    calls = {}
    with mp.Pool(max(1, min(a.threads, len(jobs)))) as pool:
        for part in pool.imap_unordered(genotype_chrom, jobs):
            calls.update(part)

    haploid_xy = a.sex == "XY"
    n_called = n_het_haploid = 0
    with open(a.out, "w") as out:
        out.write(f"# Raw genotype data in 23andMe format, written by genomeportrait v2 from whole-genome\n"
                  f"# sequencing ({a.sample}, {datetime.date.today().isoformat()}). This is NOT a 23andMe file.\n"
                  f"# Sites: Illumina Global Screening Array v3.0 SNPs. Positions: GRCh37. Bases: plus strand.\n"
                  f"# '--' = no genotype (no reads, or not callable as a SNP).\n"
                  f"# rsid\tchromosome\tposition\tgenotype\n")
        for k, (name, c, p37, p38, flip) in enumerate(rows):
            g = calls.get(k, "--")
            if g != "--" and flip:
                g = "".join(sorted(COMP[b] for b in g))
            haploid = c == "MT" or (haploid_xy and (c == "Y" or (c == "X" and not any(s <= p38 <= e for s, e in PAR38))))
            if haploid and g != "--":
                if g[0] == g[1]:
                    g = g[0]
                else:
                    g = "--"
                    n_het_haploid += 1
            if g != "--":
                n_called += 1
            out.write(f"{name}\t{c}\t{p37}\t{g}\n")
    print(f"[export_23andme] {n_called:,} of {len(rows):,} sites with a genotype; "
          f"{n_het_haploid:,} heterozygous calls on haploid chromosomes written as '--'", file=sys.stderr)


if __name__ == "__main__":
    main()
