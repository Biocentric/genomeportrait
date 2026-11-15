#!/usr/bin/env bash
# Does VerifyBamID2 recover a known amount of contamination?
#
#   SAMTOOLS="singularity exec ... samtools" VERIFYBAMID2="singularity exec ... verifybamid2" \
#   PANEL=/usr/local/share/verifybamid2-2.0.3-0/resource/1000g.phase3.100k.b38.vcf.gz.dat \
#   bash tests/validate_contamination.sh <sample.cram> <contaminant.cram> <reference.fasta> <region> <workdir>
#
# Reads of a second person are mixed into a sample at known fractions; FREEMIX should track
# them. Both CRAMs should cover the same region at similar depth (e.g. chr21 of two 1000
# Genomes 30x samples). The marker panel is restricted to <region>'s chromosome, so markers
# elsewhere do not count as zero-depth sites. The panel files are read from inside the
# VerifyBamID2 container, so it must be able to see <workdir>. One chromosome holds only a
# few hundred of the 100k markers, below VerifyBamID2's own minimum, so this test runs with
# --DisableSanityCheck; the pipeline uses the genome-wide panel with the checks on.
set -euo pipefail
SAMPLE=$1; CONT=$2; FASTA=$3; REGION=$4; WORK=$5
CHROM=${REGION%%:*}
mkdir -p "$WORK"; cd "$WORK"

# panel rows (bed, mu, UD) are per marker in the same order; V is per reference sample
$VERIFYBAMID2 sh -c "cat ${PANEL}.bed" > full.bed
$VERIFYBAMID2 sh -c "cat ${PANEL}.mu"  > full.mu
$VERIFYBAMID2 sh -c "cat ${PANEL}.UD"  > full.UD
$VERIFYBAMID2 sh -c "cat ${PANEL}.V"   > panel.V
awk -v c="$CHROM" '{ print ($1 == c) }' full.bed > keep.mask
for x in bed mu UD; do
    awk 'NR == FNR { k[FNR] = $1; next } k[FNR]' keep.mask full.$x > panel.$x
done
echo "markers on $CHROM: $(wc -l < panel.bed)" >&2

printf "target_fraction\tFREEMIX\tSNPS\tAVG_DP\n"
for f in 0 0.02 0.05 0.10; do
    if [ "$f" = 0 ]; then
        mix=$SAMPLE
    else
        # fraction of contaminant reads in the mix = s / (1 + s) at equal depth
        s=$(python3 -c "print(round($f / (1 - $f), 4))")
        $SAMTOOLS view -@ 4 -T "$FASTA" -s "7${s#0}" -b -o cont.$f.bam "$CONT" "$REGION"
        $SAMTOOLS view -@ 4 -T "$FASTA" -b -o samp.bam "$SAMPLE" "$REGION"
        $SAMTOOLS merge -@ 4 -f -c -p -o merged.$f.bam samp.bam cont.$f.bam
        # one read group: VerifyBamID2 refuses files that name two samples
        $SAMTOOLS view -h merged.$f.bam | awk '!/^@RG/' \
            | $SAMTOOLS addreplacerg -r '@RG\tID:mix\tSM:mix' -m overwrite_all -O bam -o mix.$f.bam -
        $SAMTOOLS index mix.$f.bam
        mix=mix.$f.bam
    fi
    $VERIFYBAMID2 verifybamid2 --SVDPrefix panel --Reference "$FASTA" --BamFile "$mix" \
        --DisableSanityCheck --NumThread 8 --Output out.$f > vb.$f.log 2>&1
    [ -s out.$f.selfSM ] || { echo "VerifyBamID2 failed for fraction $f, see vb.$f.log" >&2; exit 1; }
    awk -v f="$f" 'NR == 2 { printf "%s\t%s\t%s\t%s\n", f, $7, $4, $6 }' out.$f.selfSM
done
