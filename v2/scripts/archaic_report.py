#!/usr/bin/env python3
"""Neanderthal and Denisovan segments from hmmix decode, summarised for the report.

    python3 archaic_report.py --sample S --decoded S.decoded.diploid.txt \\
        --mask hg38_strick_callability_mask.bed --outdir S.archaic \\
        --segments-out S.archaic_segments.raw.tsv [--reference assets/archaic/1000g_reference.tsv]

hmmix calls a stretch archaic when it holds more variants absent from Africans than the
genome's own background rate allows. In an unphased genome a segment means at least one of
the two copies is archaic there. Segments with a mean posterior of at least 0.8 are kept.

Each segment is labelled by the archaic genomes it shares variants with (as in Skov's 1000
Genomes release): Neanderthal when the best-matching Neanderthal (Altai, Vindija or
Chagyrskaya) shares more of them than the Denisovan, Denisovan when the Denisovan shares
more, 'both' when they tie, and 'unassigned' when none are shared (a segment can be
archaic by density alone, from a lineage that none of the four sequenced genomes captures).
"""
import argparse
import os
import sys

MIN_PROB = 0.8
NEANDERTALS = ("AltaiNeandertal", "Vindija33.19", "Chagyrskaya-Phalanx")
DENISOVAN = "Denisova"
SHORT = {"AltaiNeandertal": "Altai", "Vindija33.19": "Vindija", "Chagyrskaya-Phalanx": "Chagyrskaya",
         "Denisova": "Denisova"}
COLOURS = {"Neanderthal": "#3b6ea5", "Denisovan": "#d1495b", "both": "#a56a00", "unassigned": "#9a978f"}
CHROM_LEN = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973,
             145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718,
             101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468]


def label(n, d):
    if n == d == 0:
        return "unassigned"
    if n == d:
        return "both"
    return "Neanderthal" if n > d else "Denisovan"


def read_segments(path):
    segs = []
    with open(path) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        for l in fh:
            r = dict(zip(h, l.rstrip("\n").split("\t")))
            if r["state"] != "Archaic" or float(r["mean_prob"]) < MIN_PROB:
                continue
            nea = {k: int(r.get(k, 0)) for k in NEANDERTALS}
            den = int(r.get(DENISOVAN, 0))
            segs.append({
                "chrom": r["chrom"], "start": int(r["start"]), "end": int(r["end"]),
                "length": int(r["length"]), "mean_prob": float(r["mean_prob"]),
                "snps": int(r["snps"]), "archaic_shared": int(r.get("admixpopvariants", 0)),
                **nea, DENISOVAN: den, "type": label(max(nea.values()), den),
            })
    return segs


def callable_bp(mask):
    auto = {f"chr{c}" for c in range(1, 23)}
    n = 0
    with open(mask) as fh:
        for l in fh:
            f = l.split()
            if f[0] in auto:
                n += int(f[2]) - int(f[1])
    return n


def read_reference(path):
    """Per-population archaic Mb (diploid, at least one copy) from the published 1000 Genomes calls."""
    with open(path) as fh:
        rows = [l.rstrip("\n").split("\t") for l in fh if l.strip() and not l.startswith("#")]
    return [dict(zip(rows[0], r)) for r in rows[1:]]


def plot(segs, path, sample):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    fig, ax = plt.subplots(figsize=(8.6, 5.4), dpi=130)
    for i, L in enumerate(CHROM_LEN):
        y = 22 - i
        ax.barh(y, L / 1e6, height=0.55, color="#e6e3dc", edgecolor="none")
    for s in segs:
        y = 22 - (int(s["chrom"].replace("chr", "")) - 1)
        ax.barh(y, max(s["length"] / 1e6, 0.25), left=s["start"] / 1e6, height=0.55,
                color=COLOURS[s["type"]], edgecolor="none")
    ax.set_yticks(range(1, 23), [str(23 - i) for i in range(1, 23)], fontsize=8)
    ax.set_xlabel("position (Mb)", fontsize=9)
    ax.set_ylabel("chromosome", fontsize=9)
    ax.tick_params(axis="x", labelsize=8)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    present = [t for t in COLOURS if any(s["type"] == t for s in segs)]
    ax.legend(handles=[Patch(color=COLOURS[t], label=t) for t in present], frameon=False,
              fontsize=8, loc="lower right")
    ax.set_title(f"Archaic segments in {sample} (at least one copy)", fontsize=10, loc="left")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--decoded", required=True)
    ap.add_argument("--mask", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--segments-out", required=True)
    ap.add_argument("--reference", help="per-population table from the published 1000G calls")
    a = ap.parse_args()

    segs = read_segments(a.decoded)
    os.makedirs(a.outdir, exist_ok=True)
    cols = ["chrom", "start", "end", "length", "mean_prob", "snps", "archaic_shared",
            *NEANDERTALS, DENISOVAN, "type"]
    with open(a.segments_out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for s in segs:
            fh.write("\t".join(str(s[c]) for c in cols) + "\n")

    total = sum(s["length"] for s in segs)
    by_type = {t: sum(s["length"] for s in segs if s["type"] == t) for t in COLOURS}
    call = callable_bp(a.mask)
    best = {k: sum(s[k] for s in segs) for k in (*NEANDERTALS, DENISOVAN)}
    closest = max(NEANDERTALS, key=lambda k: best[k])
    names = {"AltaiNeandertal": "Altai (Siberia, ~120,000 years ago)",
             "Vindija33.19": "Vindija (Croatia, ~50,000 years ago)",
             "Chagyrskaya-Phalanx": "Chagyrskaya (Siberia, ~80,000 years ago)"}

    rows = [
        ("Archaic sequence", f"{total / 1e6:.0f} Mb in {len(segs)} segments",
         f"{100 * total / call:.1f}% of the {call / 1e9:.2f} Gb of the autosomes that can be "
         "assessed, with at least one of the two copies archaic"),
        ("Neanderthal-like", f"{by_type['Neanderthal'] / 1e6:.0f} Mb",
         "segments sharing more variants with a Neanderthal than with the Denisovan"),
        ("Denisovan-like", f"{by_type['Denisovan'] / 1e6:.1f} Mb",
         "segments sharing more variants with the Denisovan; small amounts are expected outside "
         "Oceania and Asia, because part of the Neanderthal and Denisovan genomes is shared"),
        ("Unassigned or both", f"{(by_type['unassigned'] + by_type['both']) / 1e6:.0f} Mb",
         "archaic by variant density, without variants that tell the two groups apart"),
        ("Closest Neanderthal", names[closest] if best[closest] else "none",
         "archaic variants shared, summed over segments: "
         + ", ".join(f"{SHORT[k]} {best[k]:,}" for k in (*NEANDERTALS, DENISOVAN))),
        ("Longest segment", f"{max((s['length'] for s in segs), default=0) / 1e3:.0f} kb",
         "older admixture leaves shorter pieces: recombination cuts them every generation"),
    ]
    if a.reference:
        ref = read_reference(a.reference)
        mb = total / 1e6
        near = sorted(ref, key=lambda r: abs(float(r["median_Mb"]) - mb))[:3]
        rows.append(("Typical for", ", ".join(f"{r['population']} ({float(r['median_Mb']):.0f} Mb)" for r in near),
                     "1000 Genomes populations whose median archaic sequence is closest to this "
                     "genome's (published hmmix calls, same 0.8 cut-off). Medians run from ~72 Mb "
                     "(Puerto Rico) through ~80 (Europe) and ~95 (South Asia) to ~105 (East Asia); "
                     "African genomes, the method's outgroup, carry far less"))

    with open(os.path.join(a.outdir, "summary.tsv"), "w") as fh:
        fh.write(f"# Neanderthal and Denisovan ancestry for {a.sample} (hmmix)\n")
        fh.write("measure\tvalue\tnote\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")

    top = sorted(segs, key=lambda s: -s["length"])[:10]
    with open(os.path.join(a.outdir, "longest_segments.tsv"), "w") as fh:
        fh.write(f"# the ten longest archaic segments in {a.sample}\n")
        fh.write("region\tlength_kb\ttype\tarchaic_variants\tAltai\tVindija\tChagyrskaya\tDenisova\tposterior\n")
        for s in top:
            fh.write(f"{s['chrom']}:{s['start'] / 1e6:.2f}-{s['end'] / 1e6:.2f} Mb\t{s['length'] / 1e3:.0f}\t"
                     f"{s['type']}\t{s['archaic_shared']}\t{s['AltaiNeandertal']}\t{s['Vindija33.19']}\t"
                     f"{s['Chagyrskaya-Phalanx']}\t{s['Denisova']}\t{s['mean_prob']:.2f}\n")
    if segs:
        plot(segs, os.path.join(a.outdir, "archaic_map.png"), a.sample)
    print(f"[archaic_report] {a.sample}: {total / 1e6:.0f} Mb archaic in {len(segs)} segments "
          f"({by_type['Neanderthal'] / 1e6:.0f} Mb Neanderthal-like, {by_type['Denisovan'] / 1e6:.1f} Mb Denisovan-like)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
