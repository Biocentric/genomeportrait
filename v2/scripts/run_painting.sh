#!/usr/bin/env bash
# Chromosome painting of one genome: phase with Beagle, paint with RFMix.
#
#   bash run_painting.sh <genotypes.vcf(.gz)> <reference dir> <threads> <output prefix>
#
# <genotypes>: the person's genotypes at <reference dir>/sites.pvar, hom-ref included
# (pgs_engine.py target on the gVCF writes them; chromosomes named 1..22 or chr1..chr22).
# <reference dir>: as built by scripts/paint_reference.sh. Writes <prefix>.chrN.msp.tsv
# (RFMix segments: both haplotypes' group along the chromosome) and <prefix>.chrN.rfmix.Q
# (per-haplotype group fractions on that chromosome).
#
# Beagle phases each autosome against the 1000 Genomes panel (missing genotypes are filled
# in the process); RFMix then assigns stretches of each phased haplotype to one of the four
# reference groups. Chromosomes run in parallel, a few threads each.
set -euo pipefail
GT=$(readlink -f "$1"); REF=$(readlink -f "$2"); THREADS=$3; OUT=$(readlink -f "$4")
WORK=$(readlink -f "$(mktemp -d "${OUT}.paint.XXXX")")
PAR=$(( THREADS >= 8 ? 4 : (THREADS >= 4 ? 2 : 1) ))
PER=$(( THREADS / PAR > 0 ? THREADS / PAR : 1 ))
# RFMix otherwise searches for its CRF weight per chromosome by simulating ~500 admixed
# genomes, which is 90% of its run time (20 min per chromosome on 4 threads). On this
# reference the search picked 16-23 and its score levels off from ~12; 15 is fixed instead.
CRF_WEIGHT=15

for c in $(seq 1 22); do echo "$c chr$c"; done > "$WORK/chr_names.txt"
bcftools annotate --rename-chrs "$WORK/chr_names.txt" -Oz -o "$WORK/target.vcf.gz" "$GT"
bcftools index -f -t "$WORK/target.vcf.gz"

one() {
    c=$1
    cd "$WORK"
    bcftools view -Oz -o t.chr$c.vcf.gz target.vcf.gz chr$c
    beagle gt=t.chr$c.vcf.gz ref="$REF/phase.chr$c.vcf.gz" map="$REF/plink.chr$c.map" \
        out=phased.chr$c chrom=chr$c impute=false nthreads="$PER" > beagle.chr$c.log 2>&1
    bcftools index -f -t phased.chr$c.vcf.gz
    rfmix -f phased.chr$c.vcf.gz -r "$REF/rfmix.chr$c.vcf.gz" -m "$REF/rfmix_samples.tsv" \
        -g "$REF/rfmix.chr$c.map" -o "$OUT.chr$c" --chromosome=chr$c --n-threads="$PER" \
        -w "$CRF_WEIGHT" --random-seed=1 \
        > rfmix.chr$c.log 2>&1
    rm -f "$OUT.chr$c.fb.tsv" "$OUT.chr$c.sis.tsv"   # per-site posteriors: large, not used
}
export -f one; export WORK REF OUT PER CRF_WEIGHT
seq 1 22 | xargs -P "$PAR" -I{} bash -c 'one {} || { echo "chr{} failed:" >&2; tail -5 "$WORK"/*.chr{}.log >&2; exit 255; }'
for c in $(seq 1 22); do test -s "$OUT.chr$c.msp.tsv"; done
rm -rf "$WORK"
