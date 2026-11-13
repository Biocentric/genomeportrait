#!/usr/bin/env bash
# Run the v2 polygenic-score engine outside Nextflow, step by step.
#
# This is the same sequence the VCF_POLYGENIC subworkflow runs, written out so that
# each step can be inspected, timed and compared against a pgsc_calc run.
#
#   bash tests/pgs_standalone.sh <workdir> <scorefile_dir> <panel_dir> <gvcf|target.vcf.gz> <sample>
#
#   workdir        scratch + outputs (local disk, not a FUSE mount)
#   scorefile_dir  *_hmPOS_GRCh38.txt.gz from scripts/pgs_fetch.py
#   panel_dir      extracted pgsc_HGDP+1kGP_v1 (GRCh38_HGDP+1kGP_ALL.{pgen,pvar.zst,psam} + king cutoff ids)
#   gvcf           DeepVariant gVCF; or, to reuse an existing panel-space VCF, a file ending in target.vcf.gz
#
# Environment: PLINK2_SIF (plink2 Singularity image), THREADS (default 16), MIN_GQ (default 0),
#              BINDS (extra -B paths, comma separated), TARGET_OPTS (extra options for `target`).
set -euo pipefail

WORK=$1; SF=$2; PANEL=$3; INPUT=$4; SAMPLE=$5
HERE=$(cd "$(dirname "$0")/.." && pwd)
ENGINE="python3 $HERE/scripts/pgs_engine.py"
THREADS=${THREADS:-16}
MIN_GQ=${MIN_GQ:-0}
PLINK2_SIF=${PLINK2_SIF:?set PLINK2_SIF to a plink2 Singularity image}
BINDS=${BINDS:-}

mkdir -p "$WORK"/{ref,scores,ref_scores,target,target_scores}
cd "$WORK"     # singularity cannot start from a FUSE cwd
bind="-B $WORK -B $PANEL -B $SF"; for b in ${BINDS//,/ }; do bind="$bind -B $b"; done
plink2() { singularity exec $bind "$PLINK2_SIF" plink2 --threads "$THREADS" "$@"; }
stamp() { echo "[$(date +%T)] $*"; }

PFX="$PANEL/GRCh38_HGDP+1kGP_ALL"
KING="$PANEL/GRCh38_HGDP+1kGP.king.cutoff.out.id"

# ---- sample-independent: done once per panel -----------------------------------------
if [ ! -s sites.bed1 ]; then
    stamp "1. positions scored by the panel"
    $ENGINE sites --scorefiles "$SF" --out sites.bed1 --threads "$THREADS"
fi

if [ ! -s ref/ref.pgen ]; then
    stamp "2. subset the reference panel to those positions, unrelated samples only"
    plink2 --pfile "$PFX" vzs --extract bed1 sites.bed1 --remove "$KING" --autosome \
           --make-pgen --out ref/ref
    plink2 --pfile ref/ref --freq --out ref/ref
fi

if [ ! -s ref/ref_pcs.sscore ]; then
    stamp "3. reference PCA on common, LD-pruned SNPs"
    plink2 --pfile ref/ref --maf 0.05 --geno 0.01 --snps-only just-acgt --max-alleles 2 \
           --exclude bed1 "$HERE/assets/pgs/pca_exclude_grch38.bed1" \
           --indep-pairwise 1000kb 0.1 --out ref/prune
    plink2 --pfile ref/ref --extract ref/prune.prune.in --freq counts --pca 10 allele-wts --out ref/pca
    # Project the reference with the same command later used for the sample, so both
    # sit on one scale (plink2's .eigenvec is scaled differently from a --score projection).
    A1=$(awk 'NR==1{for(i=1;i<=NF;i++) if($i=="A1") print i}' ref/pca.eigenvec.allele)
    P1=$(awk 'NR==1{for(i=1;i<=NF;i++) if($i=="PC1") print i}' ref/pca.eigenvec.allele)
    echo "$A1 $P1 $((P1 + 9))" > ref/pca.cols
    plink2 --pfile ref/ref --read-freq ref/pca.acount \
           --score ref/pca.eigenvec.allele 2 "$A1" header-read no-mean-imputation variance-standardize \
           --score-col-nums "$P1-$((P1 + 9))" --out ref/ref_pcs
fi

if [ ! -s match_stats.tsv ]; then
    stamp "4. match every score to the reference by chr:pos:ref:alt"
    $ENGINE match --scorefiles "$SF" --pvar ref/ref.pvar --afreq ref/ref.afreq --outdir scores --stats match_stats.tsv --threads "$THREADS"
fi

if [ ! -e ref_scores/.done ]; then
    stamp "5. score the reference samples"
    for f in scores/*.score; do
        id=$(basename "$f" .score)
        [ -s "ref_scores/$id.sscore" ] && continue
        ( singularity exec $bind "$PLINK2_SIF" plink2 --threads 4 --pfile ref/ref \
              --score "$f" 1 2 3 header cols=+scoresums --out "ref_scores/$id" > "ref_scores/$id.out" 2>&1 \
          || echo "reference scoring failed: $id" >&2 ) &
        while [ "$(jobs -r | wc -l)" -ge $((THREADS / 4)) ]; do wait -n; done
    done
    wait
    touch ref_scores/.done
fi

# ---- per sample --------------------------------------------------------------------------
if [ ! -s target/target.pgen ]; then
    case "$INPUT" in
        *target.vcf.gz) stamp "6. reusing panel-space VCF $INPUT"; TVCF=$INPUT; bind="$bind -B $(dirname "$INPUT")" ;;
        *)  stamp "6. rebuild the sample in the reference's allele space (min GQ $MIN_GQ)"
            TVCF=target/$SAMPLE.target.vcf.gz
            $ENGINE target --gvcf "$INPUT" --pvar ref/ref.pvar --sample "$SAMPLE" --min-gq "$MIN_GQ" ${TARGET_OPTS:-} \
                    --out "$TVCF" --stats target/target_stats.tsv ;;
    esac
    plink2 --vcf "$TVCF" --autosome --extract ref/ref.pvar --make-pgen --out target/target
fi

stamp "7. project the sample onto the reference PCs"
read -r A1 P1 P2 < ref/pca.cols
plink2 --pfile target/target --read-freq ref/pca.acount \
       --score ref/pca.eigenvec.allele 2 "$A1" header-read no-mean-imputation variance-standardize \
       --score-col-nums "$P1-$P2" --out target/target_pcs

stamp "8. score the sample (missing genotypes take the reference mean)"
for f in scores/*.score; do
    id=$(basename "$f" .score)
    ( singularity exec $bind "$PLINK2_SIF" plink2 --threads 2 --pfile target/target --read-freq ref/ref.afreq \
          --score "$f" 1 2 3 header cols=+scoresums --out "target_scores/$id" > "target_scores/$id.out" 2>&1 \
      || echo "sample scoring failed: $id" >&2 ) &
    while [ "$(jobs -r | wc -l)" -ge $((THREADS / 2)) ]; do wait -n; done
done
wait

stamp "9. ancestry-adjusted Z and percentiles"
$ENGINE normalise --ref-scores ref_scores --target-scores target_scores \
        --ref-pcs ref/ref_pcs.sscore --target-pcs target/target_pcs.sscore \
        --psam ref/ref.psam --match-stats match_stats.tsv \
        --out "$SAMPLE.pgs_scores.tsv" --pop-out "$SAMPLE.pgs_population.tsv"
stamp "done: $WORK/$SAMPLE.pgs_scores.tsv"
