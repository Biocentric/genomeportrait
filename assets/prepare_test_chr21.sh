#!/usr/bin/env bash
#
# prepare_test_chr21.sh — fetch real 30x WGS data and slice chromosome 21 to make a
# self-contained test dataset for the `test_30x` profile.
#
# Produces (under ./testdata):
#   NA12878.chr21.30x_R1.fastq.gz / _R2.fastq.gz   (paired reads → tests alignment)
#   NA12878.chr21.30x.bam                          (also handy for a BAM-start run)
#   samplesheet_chr21_30x.csv                      (consumed by -profile test_30x)
#
# Requirements: samtools (>=1.15) + curl/wget on PATH. ~a few GB of network + disk.
#
# The default source is the NYGC/1000G 2504 high-coverage CRAM for NA12878. CRAM slicing
# streams only chr21, so you do NOT download the whole genome. If the URL has moved,
# pass a replacement CRAM (any GRCh38 30x WGS sample) as the first argument:
#   bash assets/prepare_test_chr21.sh <cram_url> [reference_fasta_url]
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${SCRIPT_DIR}/../testdata"
mkdir -p "$OUT"

CRAM_URL="${1:-https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/data/ERR3239334/NA12878.final.cram}"
# GRCh38 reference the CRAM was aligned against (needed to decode a CRAM slice).
REF_URL="${2:-https://storage.googleapis.com/genomics-public-data/references/GRCh38/v0/Homo_sapiens_assembly38.fasta}"
REGION="chr21"
THREADS="${THREADS:-4}"

echo ">> Using CRAM : $CRAM_URL"
echo ">> Region     : $REGION"
echo ">> Output dir : $OUT"

# samtools can read both the CRAM and its reference over https and slice a region remotely.
echo ">> Slicing $REGION (streaming; no full-genome download)..."
samtools view -@ "$THREADS" -T "$REF_URL" -b \
    "$CRAM_URL" "$REGION" \
    -o "$OUT/NA12878.chr21.30x.unsorted.bam"

echo ">> Sorting by name for FASTQ extraction..."
samtools sort -n -@ "$THREADS" -o "$OUT/NA12878.chr21.30x.namesort.bam" \
    "$OUT/NA12878.chr21.30x.unsorted.bam"

echo ">> Writing paired FASTQs..."
samtools fastq -@ "$THREADS" \
    -1 "$OUT/NA12878.chr21.30x_R1.fastq.gz" \
    -2 "$OUT/NA12878.chr21.30x_R2.fastq.gz" \
    -0 /dev/null -s /dev/null -n \
    "$OUT/NA12878.chr21.30x.namesort.bam"

echo ">> Writing coordinate-sorted BAM (for optional BAM-start runs)..."
samtools sort -@ "$THREADS" -o "$OUT/NA12878.chr21.30x.bam" \
    "$OUT/NA12878.chr21.30x.unsorted.bam"
samtools index "$OUT/NA12878.chr21.30x.bam"
rm -f "$OUT/NA12878.chr21.30x.unsorted.bam" "$OUT/NA12878.chr21.30x.namesort.bam"

echo ">> Writing samplesheet..."
cat > "$OUT/samplesheet_chr21_30x.csv" <<CSV
sample,fastq_1,fastq_2,bam,sex,platform
NA12878_chr21_30x,${OUT}/NA12878.chr21.30x_R1.fastq.gz,${OUT}/NA12878.chr21.30x_R2.fastq.gz,,XX,illumina
CSV

echo ">> Done. Run:  nextflow run . -profile test_30x,singularity --outdir results_test30x"
