#!/usr/bin/env python3
"""Chromosome painting from RFMix: genome-wide shares, ancestry tracts and a painted karyogram.

    python3 painting_report.py --sample S --msp-prefix S --outdir S.painting \\
        --segments-out S.painting_segments.raw.tsv

Reads <prefix>.chr1..22.msp.tsv. Each chromosome comes as two phased copies; which copy
is maternal or paternal is unknown, and copy "A" of one chromosome has nothing to do with
copy "A" of the next (each chromosome is phased on its own). Shares are therefore given for
the genome as a whole (both copies), in centimorgans, which weighs regions by how much they
recombine rather than by their length in bases.

Tracts: neighbouring windows with the same group are merged. Long tracts of a minority group
point to recent mixing (they are cut by recombination each generation); short scattered ones
in a genome that is otherwise one group are mostly noise or very old shared ancestry.
"""
import argparse
import os
import sys

LABEL = {"African": "African", "European": "European", "East_Asian": "East Asian", "South_Asian": "South Asian"}
COLOURS = {"African": "#c0803a", "European": "#3b6ea5", "East_Asian": "#1f8a70", "South_Asian": "#a14c9c"}
CHROM_LEN = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973,
             145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718,
             101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468]


def read_msp(path, sample):
    """RFMix segments for one sample (a query VCF may hold several): its two copies' groups."""
    with open(path) as fh:
        codes = fh.readline().lstrip("#").split(":", 1)[1].split()
        names = {v: k for k, v in (c.split("=") for c in codes)}
        head = fh.readline().lstrip("#").rstrip("\n").split("\t")
        i, j = head.index(f"{sample}.0"), head.index(f"{sample}.1")
        rows = []
        for l in fh:
            f = l.rstrip("\n").split("\t")
            rows.append((f[0], int(f[1]), int(f[2]), float(f[3]), float(f[4]), names[f[i]], names[f[j]]))
    return rows


def tracts(rows, copy):
    """Merge consecutive windows of one copy with the same group: (chrom, start, end, cM, group)."""
    out = []
    for chrom, s, e, sg, eg, a, b in rows:
        g = a if copy == 0 else b
        if out and out[-1][0] == chrom and out[-1][4] == g:
            out[-1][2] = e
            out[-1][3] += eg - sg
        else:
            out.append([chrom, s, e, eg - sg, g])
    return out


def plot(per_chrom, path, sample):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    fig, ax = plt.subplots(figsize=(8.6, 6.4), dpi=130)
    h = 0.34
    for i, L in enumerate(CHROM_LEN):
        y = 22 - i
        for k, off in ((0, h / 2 + 0.02), (1, -h / 2 - 0.02)):
            ax.barh(y + off, L / 1e6, height=h, color="#e6e3dc", edgecolor="none")
            for chrom, s, e, _, g in per_chrom.get(i + 1, ([], []))[k]:
                ax.barh(y + off, (e - s) / 1e6, left=s / 1e6, height=h, color=COLOURS[g], edgecolor="none")
    ax.set_yticks(range(1, 23), [str(23 - i) for i in range(1, 23)], fontsize=8)
    ax.set_xlabel("position (Mb)", fontsize=9)
    ax.set_ylabel("chromosome (two copies each)", fontsize=9)
    ax.tick_params(axis="x", labelsize=8)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    present = [g for g in COLOURS if any(t[4] == g for c in per_chrom.values() for k in c for t in k)]
    ax.legend(handles=[Patch(color=COLOURS[g], label=LABEL[g]) for g in present], frameon=False,
              fontsize=8, loc="lower right")
    ax.set_title(f"Chromosome painting of {sample}", fontsize=10, loc="left")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True, help="report label")
    ap.add_argument("--vcf-sample", help="sample name in the RFMix output, if not --sample")
    ap.add_argument("--msp-prefix", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--segments-out", required=True)
    a = ap.parse_args()

    per_chrom, all_tracts = {}, []
    cm = {g: 0.0 for g in LABEL}
    for c in range(1, 23):
        rows = read_msp(f"{a.msp_prefix}.chr{c}.msp.tsv", a.vcf_sample or a.sample)
        t = (tracts(rows, 0), tracts(rows, 1))
        per_chrom[c] = t
        for k, copy in enumerate(t):
            for chrom, s, e, length, g in copy:
                cm[g] += length
                all_tracts.append((chrom, s, e, length, g, "AB"[k]))
    total = sum(cm.values())
    share = {g: 100 * cm[g] / total for g in cm}
    main_group = max(share, key=share.get)

    os.makedirs(a.outdir, exist_ok=True)
    with open(a.segments_out, "w") as fh:
        fh.write("chrom\tstart\tend\tcM\tgroup\tcopy\n")
        for chrom, s, e, length, g, k in all_tracts:
            fh.write(f"{chrom}\t{s}\t{e}\t{length:.2f}\t{g}\t{k}\n")

    with open(os.path.join(a.outdir, "shares.tsv"), "w") as fh:
        fh.write(f"# chromosome painting of {a.sample}: share of the genome (both copies, by cM)\n")
        fh.write("group\tshare\ttracts\tlongest_tract_cM\tmedian_tract_cM\tnote\n")
        for g in sorted(LABEL, key=lambda x: -share[x]):
            mine = sorted(t[3] for t in all_tracts if t[4] == g)
            if not mine:
                fh.write(f"{LABEL[g]}\t0%\t0\t\t\t\n")
                continue
            note = ""
            if g != main_group and share[g] < 2:
                note = "below ~2% is within the noise of the method for a single genome"
            elif g != main_group and max(mine) > 50:
                note = "long tracts: mixing within the last few generations"
            fh.write(f"{LABEL[g]}\t{share[g]:.1f}%\t{len(mine)}\t{max(mine):.0f}\t{mine[len(mine) // 2]:.1f}\t{note}\n")

    with open(os.path.join(a.outdir, "about.tsv"), "w") as fh:
        fh.write("# how the painting is made\naspect\tdetail\n")
        fh.write("reference\t400 genomes from the 1000 Genomes Project: five unadmixed populations each "
                 "for African, European, East Asian and South Asian ancestry\n")
        fh.write("method\tBeagle 5.4 phases the genome into its two chromosome copies; RFMix v2 assigns "
                 "each stretch of each copy to the reference group it resembles most\n")
        fh.write("limits\tOnly these four groups exist in the reference. Ancestry from elsewhere "
                 "(Native American, Oceanian, Middle Eastern, North African...) is painted as the "
                 "nearest of the four; the copies in the picture are not maternal and paternal\n")

    plot(per_chrom, os.path.join(a.outdir, "painting.png"), a.sample)
    print(f"[painting_report] {a.sample}: " + ", ".join(f"{LABEL[g]} {share[g]:.1f}%" for g in
          sorted(LABEL, key=lambda x: -share[x])), file=sys.stderr)


if __name__ == "__main__":
    main()
