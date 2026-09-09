#!/usr/bin/env python3
"""Build an ExpansionHunter variant catalog for the long CODIS/ESS STR loci.

HipSTR cannot genotype loci whose reference tract exceeds what a 2x150 bp read pair can
span (~71 bp here) and silently omits them rather than emitting a filtered record.
ExpansionHunter combines spanning, flanking and in-repeat reads against a sequence graph,
so it is not capped by read length.

These loci are COMPOUND repeats (vWA is (AGAT)n(AGAC)n(AGAT)n, D2S1338 is (GGAA)n(GGCA)n,
...), so a single-motif structure would mis-align half the tract. Segment the reference
tract into maximal runs of its repeat unit and emit a compound LocusStructure with one
ReferenceRegion per run. The structure is verified to reconstruct the reference exactly.

Coordinates are 0-based half-open (BED style), matching both the HipSTR reference BED the
pipeline already uses and ExpansionHunter's ReferenceRegion convention.
"""
import argparse
import json
import sys
import urllib.request

MIN_COPIES = 3   # a run shorter than this is emitted as literal sequence, not a (M)* segment


def fetch_ucsc(genome, chrom, start, end):
    url = (f"https://api.genome.ucsc.edu/getData/sequence?genome={genome}"
           f";chrom={chrom};start={start};end={end}")
    with urllib.request.urlopen(url, timeout=120) as r:
        return json.load(r)["dna"].upper()


def run_at(seq, i, period):
    """Motif and copy count of the maximal exact repeat starting at position i."""
    n = len(seq)
    if i + period > n:
        return None, 0
    motif = seq[i:i + period]
    k = 0
    while i + (k + 1) * period <= n and seq[i + k * period:i + (k + 1) * period] == motif:
        k += 1
    return motif, k


def segment(seq, period):
    """Split into [(kind, text, length)]; kind is 'repeat' (text=motif) or 'literal'.

    Greedy scanning locks onto whatever phase position i lands in, which splits a long tract
    into a short literal plus a shifted motif (vWA came out as AGATG(GATA)* instead of
    (AGAT)*). Look ahead across the possible phases and take the one covering the most bases,
    so the starred segments stay as long as possible - they are what absorbs allele-length
    variation in the graph.
    """
    segs, i, n = [], 0, len(seq)
    while i < n:
        best = None
        for off in range(period):
            motif, k = run_at(seq, i + off, period)
            if k >= MIN_COPIES:
                gain = k * period - off
                if best is None or gain > best[3]:
                    best = (motif, off, k, gain)
        if best:
            motif, off, k, _ = best
            if off:                                   # bases skipped to reach the better phase
                if segs and segs[-1][0] == "literal":
                    segs[-1] = ("literal", segs[-1][1] + seq[i:i + off], segs[-1][2] + off)
                else:
                    segs.append(("literal", seq[i:i + off], off))
            segs.append(("repeat", motif, k * period))
            i += off + k * period
        else:
            if segs and segs[-1][0] == "literal":
                segs[-1] = ("literal", segs[-1][1] + seq[i], segs[-1][2] + 1)
            else:
                segs.append(("literal", seq[i], 1))
            i += 1
    return segs


def build(locus, chrom, start, end, period, genome, seq=None):
    seq = seq or fetch_ucsc(genome, chrom, start, end)
    segs = segment(seq, period)
    rebuilt = "".join(t * (l // len(t)) if k == "repeat" else t for k, t, l in segs)
    if rebuilt != seq:
        raise SystemExit(f"{locus}: structure does not reconstruct the reference tract")

    structure, regions, ids, types, pos = "", [], [], [], start
    for kind, text, length in segs:
        if kind == "repeat":
            structure += f"({text})*"
            regions.append(f"{chrom}:{pos}-{pos + length}")
            ids.append(f"{locus}_{text}_{len(regions)}")
            types.append("Repeat")
        else:
            structure += text
        pos += length

    entry = {"LocusId": locus, "LocusStructure": structure,
             "ReferenceRegion": regions, "VariantId": ids, "VariantType": types}
    if len(regions) == 1:                      # EH expects scalars for single-variant loci
        entry.update(ReferenceRegion=regions[0], VariantId=ids[0], VariantType=types[0])
    return entry, seq


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--loci", required=True,
                    help="locus:chrom:start:end:period, comma-separated (0-based half-open)")
    ap.add_argument("--genome", default="hg38")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    catalog = []
    for spec in a.loci.split(","):
        locus, chrom, start, end, period = spec.split(":")
        entry, seq = build(locus, chrom, int(start), int(end), int(period), a.genome)
        nrep = len(entry["VariantType"]) if isinstance(entry["VariantType"], list) else 1
        print(f"{locus:9s} {len(seq):3d} bp  {nrep} repeat segment(s)  {entry['LocusStructure']}",
              file=sys.stderr)
        catalog.append(entry)

    with open(a.out, "w", encoding="utf-8") as fh:
        json.dump(catalog, fh, indent=2)
        fh.write("\n")
    print(f"wrote {len(catalog)} loci to {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
