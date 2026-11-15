#!/usr/bin/env bash
# Blood-group predictions against known population frequencies.
#
#   PLINK2="singularity exec ... plink2" bash tests/validate_blood_groups.sh <HGDP+1kGP panel dir> <workdir> [samples.tsv]
#
# 1. ABO, Kell, Duffy, Kidd: the defining variants of all unrelated panel genomes go through
#    scripts/blood_groups.py; phenotype frequencies per super-population are printed, to be
#    read against published ones (e.g. Reid & Lomas-Francis, The Blood Group Antigen FactsBook).
# 2. RhD (optional, needs network): for each line of samples.tsv (sample, superpop, population,
#    CRAM URL of a 1000 Genomes 30x alignment), reads in RHD and the control windows are
#    counted exactly as RHD_DEPTH does, giving RHD copy numbers per population.
#    SAMTOOLS and FASTA must be set for this part.
set -euo pipefail
PANEL=$1; WORK=$2; SAMPLES=${3:-}
HERE=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$WORK"; cd "$WORK"
P=$PANEL/GRCh38_HGDP+1kGP_ALL
grep -v '^#' "$HERE/assets/blood_groups/sites.pvar" | cut -f3 > ids.txt
$PLINK2 --pfile "$P" vzs --remove "$PANEL/GRCh38_HGDP+1kGP.king.cutoff.out.id" --extract ids.txt \
    --export vcf bgz --out panel_bg > plink.log 2>&1
python3 "$HERE/scripts/blood_groups.py" --vcf panel_bg.vcf.gz --psam "$P.psam" --summary-out by_population.tsv
column -t -s $'\t' by_population.tsv

[ -n "$SAMPLES" ] || exit 0
echo
printf "sample\tsuperpop\tRHD_copies\tratio\n"
while IFS=$'\t' read -r s sp _pop url; do
    : > "$s.counts"
    while IFS=$'\t' read -r c a b name; do
        n=$($SAMTOOLS view -c -q 20 -F 0xD04 -T "$FASTA" "$url" "$c:$((a + 1))-$b")
        printf "%s\t%s\t%s\n" "$name" "$n" "$((b - a))" >> "$s.counts"
    done < "$HERE/assets/blood_groups/rhd_regions.bed"
    awk -F'\t' -v s="$s" -v sp="$sp" '$1 == "RHD" { r = $2 / $3 } $1 ~ /^ctrl/ { c += $2; l += $3 }
        END { x = r / (c / l); printf "%s\t%s\t%d\t%.2f\n", s, sp, (x < 0.2 ? 0 : x < 0.75 ? 1 : 2), x }' "$s.counts"
done < "$SAMPLES"
