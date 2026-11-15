#!/usr/bin/env bash
# Archaic segments from the pipeline's hmmix step on 1000 Genomes genomes, compared with the
# published hmmix calls for the same people (Skov, Zenodo 14136628, CC-BY 4.0).
#
#   PLINK2="singularity exec ... plink2" HMMIX_SIF=hmmix.sif \
#   bash tests/validate_archaic.sh <hmmix resources dir> <published segments> \
#                                  <panel pfile prefix> <samples.tsv> <workdir>
#
# <hmmix resources dir>: as stored by modules/local/hmmix_resources.nf.
# samples.tsv: sample, population (each in the panel and in the published set).
#
# Each genome is cut from the HGDP+1kGP panel used for ancestry and PGS (all SNVs, rare
# ones included) into a single-sample VCF of the sites where that person carries an ALT
# allele, which is what the pipeline hands hmmix, and run through scripts/run_hmmix.sh and
# scripts/archaic_report.py exactly as modules/local/archaic_hmmix.nf runs them.
# The published calls are per phased haplotype; their union is what a diploid run reports
# (at least one archaic copy). Per genome: archaic Mb in both, the base-pair overlap, the
# share of the published segments recovered, and agreement of the Neanderthal / Denisovan
# labels on the overlapping sequence.
set -euo pipefail
RES=$(readlink -f "$1"); PUB=$(readlink -f "$2"); PANEL=$3; SAMPLES=$(readlink -f "$4"); WORK=$5
HERE=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$WORK"; cd "$WORK"
H="singularity exec -B $WORK -B $RES -B $HERE $HMMIX_SIF"
stamp() { echo "[$(date +%T)] $*" >&2; }
blank_info() {
    gzip -dc "$1" | awk 'BEGIN { OFS = "\t" } /^#/ { print; next } { $8 = "."; print }' \
        | $H bcftools view -Oz -o "$2"
}

cut -f1 "$SAMPLES" > keep.ids
mkdir -p vcf
for c in $(seq 1 22); do
    [ -s vcf/chr$c.vcf.gz ] && continue
    $PLINK2 --pfile "$PANEL" vzs --keep keep.ids --chr $c --snps-only just-acgt --max-alleles 2 \
        --output-chr chr26 --export vcf bgz --out vcf/tmp$c > /dev/null
    # the panel's INFO holds undeclared tags, which bcftools will neither subset nor strip
    blank_info vcf/tmp$c.vcf.gz vcf/chr$c.part.vcf.gz
    mv vcf/chr$c.part.vcf.gz vcf/chr$c.vcf.gz; rm -f vcf/tmp$c.*
done
[ -s all.vcf.gz.tbi ] || { for c in $(seq 1 22); do echo vcf/chr$c.vcf.gz; done > vcf.list
    $H bcftools concat -f vcf.list -Oz -o all.vcf.gz -W=tbi 2> /dev/null; }

while read -r -u 3 s; do
    [ -s $s.archaic/summary.tsv ] && continue
    stamp "$s"
    $H bcftools view -s $s -c1 -Oz -o $s.vcf.gz -W=tbi all.vcf.gz
    $H bash "$HERE/scripts/run_hmmix.sh" $s.vcf.gz "$RES" 8 $s > $s.log 2>&1
    $H python3 "$HERE/scripts/archaic_report.py" --sample $s --decoded $s.decoded.diploid.txt \
        --mask "$RES/hg38_strick_callability_mask.bed" --outdir $s.archaic \
        --segments-out $s.archaic_segments.raw.tsv --reference "$HERE/assets/archaic/1000g_reference.tsv" >> $s.log 2>&1
done 3< keep.ids

python3 - "$SAMPLES" "$PUB" <<'EOF'
import collections, statistics, sys
samples = [l.split()[0:2] for l in open(sys.argv[1]) if l.strip()]
want = {s for s, _ in samples}

def union(iv):
    out = []
    for c, a, b, *rest in sorted(iv):
        if out and out[-1][0] == c and a <= out[-1][2]:
            out[-1][2] = max(out[-1][2], b)
        else:
            out.append([c, a, b])
    return out

def overlap(x, y):
    by = collections.defaultdict(list)
    for c, a, b, *_ in y:
        by[c].append((a, b))
    return sum(max(0, min(b, d) - max(a, e)) for c, a, b, *_ in x for e, d in by[c])

def label(n, d):
    return "none" if n == d == 0 else "Both" if n == d else "Neanderthal" if n > d else "Denisova"

pub = collections.defaultdict(list)
with open(sys.argv[2]) as fh:
    h = fh.readline().rstrip("\n").split("\t")
    for l in fh:
        r = dict(zip(h, l.rstrip("\n").split("\t")))
        if r["name"] in want and float(r["mean_prob"]) >= 0.8:
            pub[r["name"]].append((r["chrom"], int(r["start"]), int(r["end"]), r["ND_type"]))

print("sample\tpop\tpublished_Mb\tpipeline_Mb\toverlap_Mb\trecovered_%\tprecision_%\tlabel_agree_%")
tot = []
for s, pop in samples:
    mine = []
    with open(f"{s}.archaic_segments.raw.tsv") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        for l in fh:
            r = dict(zip(h, l.rstrip("\n").split("\t")))
            mine.append((r["chrom"], int(r["start"]), int(r["end"]), r["type"]))
    p, m = union(pub[s]), union(mine)
    lp = sum(b - a for _, a, b in p); lm = sum(b - a for _, a, b in m); ov = overlap(p, m)
    # label agreement, base pairs: published Neanderthal/Denisova vs ours, where both labelled
    agree = tot_lab = 0
    for c, a, b, t in pub[s]:
        if t not in ("Neanderthal", "Denisova"):
            continue
        for c2, a2, b2, t2 in mine:
            if c2 == c and t2 in ("Neanderthal", "Denisovan"):
                o = max(0, min(b, b2) - max(a, a2))
                tot_lab += o
                agree += o if t[:4] == t2[:4] else 0
    tot.append((lp, lm))
    print(f"{s}\t{pop}\t{lp/1e6:.1f}\t{lm/1e6:.1f}\t{ov/1e6:.1f}\t{100*ov/max(lp,1):.0f}\t"
          f"{100*ov/max(lm,1):.0f}\t{100*agree/max(tot_lab,1):.0f}")
if len(tot) > 2:
    print(f"correlation of archaic Mb across genomes: r = {statistics.correlation(*zip(*tot)):.2f}")
EOF
