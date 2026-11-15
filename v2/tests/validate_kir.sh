#!/usr/bin/env bash
# KIR gene content from T1K on 1000 Genomes genomes, checked for internal consistency and
# against known population frequencies (there is no openly available KIR truth set for
# these samples).
#
#   SAMTOOLS="singularity exec ... samtools" T1K_SIF=t1k.sif FASTA=GRCh38.fa \
#   bash tests/validate_kir.sh <t1k KIR index dir> <samples.tsv> <workdir>
#
# samples.tsv: sample, superpop, population, CRAM URL (1000 Genomes 30x).
# Checks: the framework genes KIR3DL3, KIR3DP1, KIR2DL4 and KIR3DL2 are found in every
# genome; an allele counts only with T1K quality > 0 (quality 0 calls are noise from
# paralogous reads); KIR2DL2/KIR2DL3 and KIR3DL1/KIR3DS1, alleles of one locus each, never exceed two
# copies together; per-population carrier frequencies of each gene are printed.
set -euo pipefail
IDX=$1; SAMPLES=$2; WORK=$3
REGION=chr19:54700000-54900000
mkdir -p "$WORK"; cd "$WORK"
SEQ=$(ls "$IDX"/*_dna_seq.fa | head -1)
B="-B $WORK -B $(dirname "$SEQ")"

type_one() {
    s=$1; url=$2
    [ -s "$s/${s}_genotype.tsv" ] && return 0
    mkdir -p "$s"
    # the chr19 ALT contigs carry KIR genes absent from the primary haplotype (KIR2DL5,
    # KIR2DS3/5, KIR3DS1...); reads placed there must be passed to T1K as well
    alts=$($SAMTOOLS view -H "$url" 2> /dev/null | sed -n 's/.*\tSN:\(chr19_[A-Za-z0-9]*_alt\).*/\1/p' | tr '\n' ' ')
    $SAMTOOLS view -b -T "$FASTA" -o "$s/kir.bam" "$url" "$REGION" $alts 2> "$s/stream.err"
    $SAMTOOLS collate -u -O "$s/kir.bam" "$s/tmp" | $SAMTOOLS fastq -1 "$s/r1.fq.gz" -2 "$s/r2.fq.gz" -0 /dev/null -s /dev/null -n - 2>/dev/null
    singularity exec $B "$T1K_SIF" run-t1k -1 "$s/r1.fq.gz" -2 "$s/r2.fq.gz" --preset kir-wgs -f "$SEQ" \
        -t 2 --od "$s" -o "$s" > "$s/t1k.log" 2>&1 || true
}
export -f type_one; export SAMTOOLS FASTA T1K_SIF SEQ B REGION
n=0
while IFS=$'\t' read -r s _sp _pop url; do
    type_one "$s" "$url" &
    n=$((n + 1)); [ $((n % 8)) -eq 0 ] && wait
done < "$SAMPLES"
wait

python3 - "$SAMPLES" <<'EOF'
import collections, os, sys
rows = [l.rstrip("\n").split("\t") for l in open(sys.argv[1])]
FRAME = ["KIR3DL3", "KIR3DP1", "KIR2DL4", "KIR3DL2"]
GENES = ["KIR2DL1", "KIR2DL2", "KIR2DL3", "KIR2DL5A", "KIR2DL5B", "KIR2DS1", "KIR2DS2", "KIR2DS3",
         "KIR2DS4", "KIR2DS5", "KIR3DL1", "KIR3DS1", "KIR2DP1"]
carriers = collections.defaultdict(collections.Counter)
n_pop = collections.Counter()
frame_missing, locus_excess, failed = [], [], []
for s, sp, pop, _url in rows:
    path = f"{s}/{s}_genotype.tsv"
    if not os.path.exists(path):
        failed.append(s); continue
    copies = {}
    for l in open(path):
        f = l.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        gene, num = f[0], int(f[1])
        # T1K: gene, number of alleles called, allele1, abundance, quality, [allele2 ...]
        alleles = [f[i] for i in (2, 5) if len(f) > i + 2 and f[i] not in (".", "") and int(f[i + 2]) > 0]
        copies[gene] = len(alleles) if num else 0
    n_pop[sp] += 1
    for g, c in copies.items():
        if c:
            carriers[sp][g] += 1
    miss = [g for g in FRAME if not copies.get(g)]
    if miss:
        frame_missing.append((s, miss))
    for a, b in (("KIR2DL2", "KIR2DL3"), ("KIR3DL1", "KIR3DS1")):
        if copies.get(a, 0) + copies.get(b, 0) > 2:
            locus_excess.append((s, a, b))
print(f"genomes typed: {sum(n_pop.values())}, failed: {failed or 'none'}")
print(f"framework genes missing: {frame_missing or 'none'}")
print(f"more than two alleles at KIR2DL2/3 or KIR3DL1/S1: {locus_excess or 'none'}")
pops = sorted(n_pop)
print("carrier frequency %   " + "  ".join(f"{p}(n={n_pop[p]})" for p in pops))
for g in FRAME + GENES:
    print(f"  {g:9s}            " + "  ".join(f"{100 * carriers[p][g] / n_pop[p]:9.0f}" for p in pops))
EOF
