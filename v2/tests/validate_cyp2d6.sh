#!/usr/bin/env bash
# CYP2D6 from Cyrius into PharmCAT, checked on 1000 Genomes samples of known structure.
#
#   CYRIUS_SIF=genomeportrait-cyrius-1.1.1.sif PHARMCAT_SIF=pharmcat-2.15.5.sif \
#   bash tests/validate_cyp2d6.sh <reference.fasta> <workdir> <sample.chr22.cram> [...]
#
# Truth (Chen et al. 2021, Pharmacogenomics J, Cyrius paper, main text): NA12878 *3/*68+*4;
# HG02465 carries *5 (deletion); HG00290 a *13+*2 tandem; NA19982 *13 and *68 on separate
# haplotypes; HG01624 a duplication; NA19239 no structural variant.
#
# Only chr22 of each sample is used (streaming whole 30x genomes is slow), so Cyrius is
# given the 62 of its 3,004 regions that lie on chr22: the CYP2D6/CYP2D7 locus and 58
# depth-normalisation regions. The pipeline runs Cyrius unmodified on the full alignment.
# Each call is then passed to PharmCAT as an outside call, to check it is accepted and
# turned into a phenotype.
set -euo pipefail
FASTA=$1; WORK=$2; shift 2
mkdir -p "$WORK"; cd "$WORK"
B="-B $(dirname "$FASTA") -B $WORK"; for d in $(for c in "$@"; do dirname "$(readlink -f "$c")"; done | sort -u); do [ "$d" = "$WORK" ] || B="$B -B $d"; done
CYR_DATA=/opt/conda/lib/python3.9/site-packages/cyrius/data/CYP2D6_region_38.bed

singularity exec $B "$CYRIUS_SIF" cat "$CYR_DATA" | awk '$1 == "chr22"' > region38_chr22.bed
for c in "$@"; do readlink -f "$c"; done > manifest.txt
singularity exec $B -B "$WORK/region38_chr22.bed:$CYR_DATA" "$CYRIUS_SIF" \
    cyrius --manifest manifest.txt --genome 38 --prefix truth --outDir cyrius --threads 6 --reference "$FASTA" \
    > cyrius.log 2>&1

printf "sample\tcyrius\tfilter\tpharmcat_cyp2d6_diplotype\tpharmcat_phenotype\n"
tail -n +2 cyrius/truth.tsv | while IFS=$'\t' read -r s g f; do
    cram=$(grep "/$s.cram\$" manifest.txt)
    printf "CYP2D6\t%s\n" "$g" > "$s.outside.tsv"
    # the PHARMCAT module's genotyping, on this chr22 slice
    singularity exec $B "$PHARMCAT_SIF" sh -c "
        bcftools query -f '%CHROM\t%POS\n' /pharmcat/pharmcat_positions.vcf | sort -k1,1 -k2,2n -u > sites.tsv
        bcftools mpileup -R sites.tsv -f $FASTA -a FORMAT/DP -Ou $cram 2>/dev/null | bcftools call -m -Oz -o $s.raw.vcf.gz
        bcftools annotate -x FORMAT/AD,FORMAT/PL -Oz -o $s.pgx.vcf.gz $s.raw.vcf.gz && bcftools index -t -f $s.pgx.vcf.gz
        python3 /pharmcat/pharmcat_vcf_preprocessor.py -vcf $s.pgx.vcf.gz -o pre_$s -bf $s
        java -jar /pharmcat/pharmcat.jar -vcf \$(ls pre_$s/$s*.preprocessed.vcf* | grep -v csi | head -1) \
            -o pc_$s -bf $s -reporterJson -po $s.outside.tsv" > pc_$s.log 2>&1 \
        || { printf "%s\t%s\t%s\tPHARMCAT FAILED\t-\n" "$s" "$g" "$f"; continue; }
    python3 - "$s" "$g" "$f" "pc_$s/$s.report.json" <<'EOF'
import json, sys
s, g, f, path = sys.argv[1:5]
rep = json.load(open(path))
dip, phen = "-", "-"
def walk(o):
    global dip, phen
    if isinstance(o, dict):
        if o.get("geneSymbol") == "CYP2D6" and "sourceDiplotypes" in o:
            d = o["sourceDiplotypes"][0] if o["sourceDiplotypes"] else {}
            dip = d.get("label", "-")
            phen = ", ".join(d.get("phenotypes", [])) or "-"
            return
        for v in o.values():
            walk(v)
    elif isinstance(o, list):
        for v in o:
            walk(v)
walk(rep)
print(f"{s}\t{g}\t{f}\t{dip}\t{phen}")
EOF
done
