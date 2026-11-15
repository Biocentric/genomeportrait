#!/usr/bin/env bash
# Archaic segments in one unphased genome with hmmix (Skov et al. 2018), as its authors run it.
#
#   bash run_hmmix.sh <sample.vcf.gz> <resources dir> <threads> <output prefix>
#
# <resources dir> is what modules/local/hmmix_resources.nf stores: the callability mask,
# the mutation-rate bed, outgroup/chrN.txt, ancestral/homo_sapiens_ancestor_N.fa and
# archaic/individuals_highcov.chrN.bcf. Writes <prefix>.decoded.diploid.txt and
# <prefix>.trained.json.
#
#   1. the variants the person carries (PASS SNVs with at least one ALT allele) inside the
#      callability mask, minus every position where the African-rich outgroup carries a
#      derived allele, polarised with the ancestral genome: create_ingroup, per autosome
#   2. a two-state HMM (human / archaic) trained on this genome's own variant density: train
#   3. posterior decoding into segments, each annotated with the variants it shares with
#      the Altai, Vindija and Chagyrskaya Neanderthals and the Denisovan: decode
#
# Only the person's own variant sites are used. hmmix's bcftools/vcftools steps, its
# per-chromosome ancestral lookup and everything else run unchanged.
set -euo pipefail
VCF=$(readlink -f "$1"); RES=$(readlink -f "$2"); THREADS=$3; OUT=$4
MASK=$RES/hg38_strick_callability_mask.bed
MUT=$RES/hg38_mutationrate.bed
WORK=$(readlink -f "$(mktemp -d "${OUT}.hmmix.XXXX")")

# grep -c, not -q: with pipefail, grep -q stopping early would SIGPIPE bcftools on a long header
CHR=$( [ "$(bcftools view -h "$VCF" | grep -c '^##contig=<ID=chr1,')" -gt 0 ] && echo chr || true)
for c in $(seq 1 22); do echo "$c chr$c"; done > "$WORK/chr_names.txt"
# INFO is dropped: hmmix subsets samples with bcftools, which refuses undefined INFO tags
bcftools view -f PASS,. -v snps -i 'GT="alt"' -Ou "$VCF" \
    | bcftools annotate -x INFO -Oz -o "$WORK/snps.vcf.gz" -W=tbi
SAMPLE=$(bcftools query -l "$WORK/snps.vcf.gz" | sed -n 1p)
# hmmix filters with vcftools 0.1.16, which stops at a VCFv4.3 header (plink2 and newer
# callers write one); the records hmmix reads are the same in 4.2
bcftools view -h "$WORK/snps.vcf.gz" | sed '1s/^##fileformat=VCFv4\.[3-9]$/##fileformat=VCFv4.2/' \
    | bcftools reheader -h /dev/stdin -o "$WORK/snps42.vcf.gz" "$WORK/snps.vcf.gz"
mv "$WORK/snps42.vcf.gz" "$WORK/snps.vcf.gz"; bcftools index -f -t "$WORK/snps.vcf.gz"
for c in $(seq 1 22); do
    bcftools view -Ou "$WORK/snps.vcf.gz" "${CHR}$c" \
        | bcftools annotate --rename-chrs "$WORK/chr_names.txt" -Oz -o "$WORK/chr$c.vcf.gz" -W=tbi
done

ingroup() {
    c=$1
    mkdir -p "$WORK/ing_$c" && cd "$WORK/ing_$c"
    hmmix create_ingroup -ind="$SAMPLE" -vcf="$WORK/chr$c.vcf.gz" -weights="$MASK" \
        -outgroup="$RES/outgroup/chr$c.txt" -ancestral="$RES/ancestral/homo_sapiens_ancestor_$c.fa" \
        -out=obs > "$WORK/ingroup.chr$c.log" 2>&1
}
export -f ingroup; export SAMPLE WORK MASK RES
seq 1 22 | xargs -P "$THREADS" -I{} bash -c 'ingroup {}'

OBS=$WORK/obs.txt
printf "chrom\tpos\tancestral_base\tgenotype\n" > "$OBS"
for c in $(seq 1 22); do tail -n +2 "$WORK/ing_$c/obs.$SAMPLE.txt"; done >> "$OBS"
n=$(( $(wc -l < "$OBS") - 1 ))
echo "[run_hmmix] $SAMPLE: $n derived variants absent from the outgroup" >&2
# a 30x genome has ~100,000-200,000 of these (fewer in African genomes, but far more than
# this); a handful means the input was not what it should be
if [ "$n" -lt 10000 ]; then
    echo "[run_hmmix] too few variants: check the VCF (single sample, GRCh38, PASS calls)" >&2
    exit 1
fi

hmmix train -obs="$OBS" -weights="$MASK" -mutrates="$MUT" -out="$OUT.trained.json" > "$WORK/train.log" 2>&1
ARCHAIC=$(for c in $(seq 1 22); do printf "%s/archaic/individuals_highcov.chr%s.bcf," "$RES" "$c"; done | sed 's/,$//')
hmmix decode -obs="$OBS" -weights="$MASK" -mutrates="$MUT" -param="$OUT.trained.json" \
    -admixpop="$ARCHAIC" -out="$OUT.decoded" > "$WORK/decode.log" 2>&1
test -s "$OUT.decoded.diploid.txt"
rm -rf "$WORK"
