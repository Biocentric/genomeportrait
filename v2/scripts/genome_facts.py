#!/usr/bin/env python3
"""Basic facts about a genome from its depth profile: coverage, sex chromosomes, mtDNA.

    python3 genome_facts.py --sample S --summary S.mosdepth.summary.txt \
        [--dist S.mosdepth.global.dist.txt] [--declared XX|XY|unknown] \
        --out S.genome_facts.tsv --sex-out S.sex.txt

Sex is inferred from the depth of chrX and chrY relative to the autosomes. Only the two
common complements are named. Anything else is reported as "undetermined" together with
the two ratios: naming a sex-chromosome aneuploidy is a clinical statement, and a depth
ratio from one short-read run is not enough to make it.

mtDNA copy number is 2 x (chrM depth / autosomal depth), i.e. mitochondrial genomes per
diploid nuclear genome. It depends strongly on the tissue sequenced.
"""
import argparse
import sys

AUTOSOMES = {f"chr{i}" for i in range(1, 23)} | {str(i) for i in range(1, 23)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--dist")
    ap.add_argument("--declared", default="unknown")
    ap.add_argument("--out", required=True)
    ap.add_argument("--sex-out", required=True)
    a = ap.parse_args()

    length, bases = {}, {}
    with open(a.summary) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0].endswith("_region"):
                continue
            length[f[0]], bases[f[0]] = int(f[1]), int(f[2])

    auto = [c for c in length if c in AUTOSOMES]
    if not auto:
        sys.exit("genome_facts: no autosomes found in the mosdepth summary")
    auto_depth = sum(bases[c] for c in auto) / sum(length[c] for c in auto)

    def depth(*names):
        for n in names:
            if n in length and length[n]:
                return bases[n] / length[n]
        return None

    x, y, m = depth("chrX", "X"), depth("chrY", "Y"), depth("chrM", "MT", "chrMT", "M")
    xr = x / auto_depth if x is not None and auto_depth else None
    yr = y / auto_depth if y is not None and auto_depth else None

    # chrX sits near 1.0 with two copies and near 0.5 with one. chrY's mean is diluted by
    # its unmappable heterochromatin, so one copy reads ~0.2-0.45 and none reads ~0-0.03
    # (stray reads in the pseudo-autosomal and X-transposed regions).
    if xr is None or yr is None:
        inferred, why = "unknown", "chrX or chrY missing from the depth summary"
    elif xr >= 0.8 and yr < 0.06:
        inferred, why = "XX", "chrX at autosomal depth, chrY essentially uncovered"
    elif 0.35 <= xr < 0.7 and 0.12 <= yr < 0.6:
        inferred, why = "XY", "chrX at half the autosomal depth, chrY covered"
    else:
        inferred, why = "unknown", "depth ratios fit neither XX nor XY cleanly"

    declared = a.declared if a.declared in ("XX", "XY") else "unknown"
    resolved = declared if declared != "unknown" else inferred
    if declared == "unknown":
        agree = "not declared in the samplesheet; inferred value used"
    elif inferred == "unknown":
        agree = "could not be checked against the data"
    elif declared == inferred:
        agree = "matches the samplesheet"
    else:
        agree = f"DISAGREES with the samplesheet ({declared}); samplesheet value used - check the sample"

    rows = [
        ("Mean autosomal depth", f"{auto_depth:.1f}x", "chr1-22, all positions"),
    ]
    if a.dist:
        cum = {}
        with open(a.dist) as fh:
            for line in fh:
                f = line.split("\t")
                if f[0] == "total":
                    cum[int(f[1])] = float(f[2])
        for d in (1, 10, 20, 30):
            if d in cum:
                rows.append((f"Genome covered at >= {d}x", f"{100 * cum[d]:.1f}%", "all contigs, including unplaced"))
    rows += [
        ("chrX depth / autosomal depth", "n/a" if xr is None else f"{xr:.2f}", "about 1.0 with two X, 0.5 with one"),
        ("chrY depth / autosomal depth", "n/a" if yr is None else f"{yr:.2f}", "about 0.2-0.45 with a Y, near 0 without"),
        ("Sex chromosomes (inferred)", "undetermined" if inferred == "unknown" else inferred, why),
        ("Sex chromosomes (used downstream)", "undetermined" if resolved == "unknown" else resolved, agree),
    ]
    if m is not None and auto_depth:
        rows.append(("mtDNA copies per cell", f"{2 * m / auto_depth:.0f}",
                     f"2 x chrM depth ({m:.0f}x) / autosomal depth; depends on the tissue sequenced"))

    with open(a.out, "w") as fh:
        fh.write(f"# genome facts for {a.sample}\n")
        fh.write("fact\tvalue\tnote\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    with open(a.sex_out, "w") as fh:
        fh.write(resolved + "\n")
    print(f"[genome_facts] depth {auto_depth:.1f}x, X ratio {xr}, Y ratio {yr}, "
          f"inferred {inferred}, declared {declared}, used {resolved}", file=sys.stderr)


if __name__ == "__main__":
    main()
