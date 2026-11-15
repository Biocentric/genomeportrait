#!/usr/bin/env bash
# Does ALT-aware alignment keep reads at ALT-covered loci on the primary assembly?
#
#   BWA="singularity exec <bwa-mem2+samtools image>" CYRIUS_SIF=<cyrius image> \
#   bash tests/validate_alt_aware.sh <GRCh38 fasta with ALTs> <bwa .alt file> <sample.chr22.cram> <workdir>
#
# chr22 and its nine ALT contigs are extracted from the full reference, with the matching
# lines of the .alt file. The sample's chr22 reads are aligned twice, with and without the
# .alt file next to the index (bwa-mem2 reads <prefix>.alt and then maps ALT-aware). The
# CYP2D6/CYP2D7 locus has an ALT haplotype (chr22_KI270879v1_alt): without ALT-aware
# mapping its reads split between the two at MAPQ 0, and Cyrius cannot call CYP2D6.
set -euo pipefail
FASTA=$1; ALT=$2; CRAM=$3; WORK=$4
mkdir -p "$WORK"; cd "$WORK"
CONTIGS="chr22 $(cut -f1 "$ALT" | grep '^chr22_' | sort -u | tr '\n' ' ')"
LOCUS=chr22:42123192-42145745
stamp() { echo "[$(date +%T)] $*" >&2; }

if [ ! -s ref.fa.bwt.2bit.64 ]; then
    stamp "reference: $CONTIGS"
    $BWA samtools faidx "$FASTA" $CONTIGS > ref.fa
    $BWA samtools faidx ref.fa
    $BWA bwa-mem2 index ref.fa > index.log 2>&1
fi
awk -v keep="$CONTIGS" 'BEGIN { n = split(keep, k, " "); for (i = 1; i <= n; i++) want[k[i]] = 1 }
     /^@/ { next } ($1 in want)' "$ALT" > alt_lines.sam

if [ ! -s r1.fq.gz ]; then
    stamp "reads from $CRAM"
    $BWA samtools collate -@ 8 -u -O --reference "$FASTA" "$CRAM" tmpc \
        | $BWA samtools fastq -@ 4 -1 r1.fq.gz -2 r2.fq.gz -0 /dev/null -s /dev/null -n - 2>/dev/null
fi

for mode in plain altaware; do
    rm -f ref.fa.alt
    [ "$mode" = altaware ] && cp alt_lines.sam ref.fa.alt
    stamp "align ($mode)"
    $BWA bash -c "bwa-mem2 mem -t 16 -R '@RG\tID:S\tSM:S' ref.fa r1.fq.gz r2.fq.gz 2> $mode.bwa.log \
        | samtools sort -@ 8 -o $mode.bam - && samtools index $mode.bam"
    grep -m1 -i "alt" $mode.bwa.log >&2 || true
done
rm -f ref.fa.alt

printf "mode\tlocus_reads\tlocus_MAPQ0\tKI270879_alt_reads\n"
for mode in plain altaware; do
    printf "%s\t%s\t%s\t%s\n" "$mode" \
        "$($BWA samtools view -c $mode.bam $LOCUS)" \
        "$($BWA samtools view -c -e 'mapq==0' $mode.bam $LOCUS)" \
        "$($BWA samtools view -c $mode.bam chr22_KI270879v1_alt)"
done

if [ -n "${CYRIUS_SIF:-}" ]; then
    DATA=/opt/conda/lib/python3.9/site-packages/cyrius/data/CYP2D6_region_38.bed
    singularity exec "$CYRIUS_SIF" cat "$DATA" | awk '$1 == "chr22"' > region38_chr22.bed
    ls "$PWD"/plain.bam "$PWD"/altaware.bam > manifest.txt
    singularity exec -B "$PWD" -B "$PWD/region38_chr22.bed:$DATA" "$CYRIUS_SIF" \
        cyrius --manifest manifest.txt --genome 38 --prefix alt --outDir cyrius --threads 4 > cyrius.log 2>&1
    echo "Cyrius (chr22 normalisation regions only):"
    cat cyrius/alt.tsv
fi
