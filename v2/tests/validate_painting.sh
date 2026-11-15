#!/usr/bin/env bash
# Chromosome painting of 1000 Genomes genomes that were kept out of both references.
#
#   PLINK2="singularity exec ... plink2" PAINT_SIF=painting.sif \
#   bash tests/validate_painting.sh <painting reference dir> <panel pfile prefix> <workdir>
#
# The hold-out genomes (assets/painting/holdout.tsv: 2 each of GBR, YRI, CHB, GIH, and 6
# each of the admixed ASW (African American) and ACB (African Caribbean)) are absent from
# Beagle's and RFMix's references. Their genotypes come from the HGDP+1kGP panel (gnomAD's
# own calls of the same genomes, unphased), and go through scripts/run_painting.sh and
# scripts/painting_report.py as the pipeline runs them.
# Expected: the unadmixed genomes painted almost entirely as their own group; ASW around
# three-quarters African and one quarter European, ACB around 85-90% African (1000 Genomes
# Consortium 2015, global ADMIXTURE estimates), with long European tracts (recent mixing).
set -euo pipefail
REF=$(readlink -f "$1"); PANEL=$2; WORK=$3
HERE=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$WORK"; cd "$WORK"
P="singularity exec -B $WORK -B $REF -B $HERE $PAINT_SIF"
stamp() { echo "[$(date +%T)] $*" >&2; }

grep -v -e '^#' -e '^sample' "$HERE/assets/painting/holdout.tsv" > holdout.tsv
cut -f1 holdout.tsv > keep.ids
grep -v '^#' "$REF/sites.pvar" | cut -f3 > sites.ids
if [ ! -s all.vcf.gz ]; then
    $PLINK2 --pfile "$PANEL" vzs --keep keep.ids --extract sites.ids --export vcf bgz --out raw > /dev/null
    # the panel's INFO holds undeclared tags, which bcftools will not subset
    gzip -dc raw.vcf.gz | awk 'BEGIN { OFS = "\t" } /^#/ { print; next } { $8 = "."; print }' \
        | $P bcftools view -Oz -o all.vcf.gz
    rm -f raw.vcf.gz
fi
echo "sites in the reference: $(wc -l < sites.ids); found in the panel: $($P bcftools view -H all.vcf.gz | wc -l)" >&2

# all hold-outs in one run: RFMix trains its random forests on the reference once per run,
# so painting 20 genomes together costs little more than one. RFMix paints each query
# haplotype on its own; in Beagle the other hold-outs add 40 haplotypes to the 6,364 of the
# reference, too few to change the phasing
mkdir -p rfmix
if [ ! -s rfmix/holdout.chr22.msp.tsv ]; then
    stamp "phasing and painting $(wc -l < keep.ids) genomes"
    $P bash "$HERE/scripts/run_painting.sh" all.vcf.gz "$REF" 32 rfmix/holdout > paint.log 2>&1
fi
while read -r -u 3 s pop; do
    $P python3 "$HERE/scripts/painting_report.py" --sample $s --msp-prefix rfmix/holdout --outdir $s.painting \
        --segments-out $s.painting_segments.raw.tsv 2>> paint.log
done 3< holdout.tsv

python3 - <<'EOF'
import collections
rows = [l.split() for l in open("holdout.tsv")]
own = {"GBR": "European", "YRI": "African", "CHB": "East Asian", "GIH": "South Asian"}
print("sample\tpop\tAfrican\tEuropean\tEast_Asian\tSouth_Asian\tlongest_minority_tract_cM")
by_pop = collections.defaultdict(list)
for s, pop in rows:
    sh, longest = {}, {}
    for l in open(f"{s}.painting/shares.tsv"):
        if l.startswith(("#", "group")):
            continue
        f = l.rstrip("\n").split("\t")
        sh[f[0]] = float(f[1].rstrip("%"))
        longest[f[0]] = float(f[3]) if f[3] else 0.0
    main = max(sh, key=sh.get)
    minor = max((longest[g] for g in sh if g != main), default=0)
    by_pop[pop].append(sh)
    print(f"{s}\t{pop}\t" + "\t".join(f"{sh.get(g, 0):.1f}" for g in ("African", "European", "East Asian", "South Asian"))
          + f"\t{minor:.0f}")
print()
for pop, lst in by_pop.items():
    mean = {g: sum(x.get(g, 0) for x in lst) / len(lst) for g in ("African", "European", "East Asian", "South Asian")}
    print(f"{pop} mean: " + ", ".join(f"{g} {v:.1f}%" for g, v in mean.items() if v >= 0.5))
EOF
