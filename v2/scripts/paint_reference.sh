#!/usr/bin/env bash
# Reference for chromosome painting, built once from the 1000 Genomes 30x phased panel.
#
#   bash paint_reference.sh <panel base URL> <samples ped URL> <genetic map zip URL> \
#        <holdout.tsv> <threads> <outdir>
#
# Per autosome (streamed, MD5-checked against the panel's manifest, never stored whole):
#   - biallelic SNPs with panel allele frequency 5-95%, thinned to at most one per 3 kb
#     (~900,000 sites genome-wide): enough for painting, small enough for RFMix
#   - phase.chrN.vcf.gz: all panel genomes except the hold-out list, for Beagle phasing
#   - rfmix.chrN.vcf.gz: the RFMix reference, 20 unrelated genomes from each of five
#     unadmixed populations per continental group (African, European, East Asian, South
#     Asian), hold-outs excluded
# plus sites.pvar (the sites, for reading the person's genotypes from the gVCF),
# rfmix_samples.tsv (sample -> group), genetic maps for Beagle (plink.chrN.map) and RFMix
# (rfmix.chrN.map). The hold-out genomes (assets/painting/holdout.tsv) are kept out of
# both references so tests/validate_painting.sh can paint them honestly.
set -euo pipefail
BASE=$1; PED=$2; MAPZIP=$3; HOLDOUT=$(readlink -f "$4"); THREADS=$5; OUT=$6
mkdir -p "$OUT"; cd "$OUT"
stamp() { echo "[paint_reference $(date +%T)] $*" >&2; }

curl -sfL --retry 5 -o manifest.txt "$BASE/20220804_manifest.txt"
curl -sfL --retry 5 -o samples.ped "$PED"
curl -sfL --retry 5 -o maps.zip "$MAPZIP"
python3 - "$HOLDOUT" <<'EOF'
import sys, zipfile
GROUPS = {"African": ["YRI", "LWK", "ESN", "MSL", "GWD"], "European": ["CEU", "GBR", "FIN", "IBS", "TSI"],
          "East_Asian": ["CHB", "CHS", "JPT", "KHV", "CDX"], "South_Asian": ["GIH", "PJL", "BEB", "STU", "ITU"]}
PER_POP = 20
hold = {l.split()[0] for l in open(sys.argv[1]) if l.strip() and not l.startswith(("#", "sample"))}
rows = [l.split() for l in open("samples.ped")][1:]
children = {r[1] for r in rows if r[2] != "0" or r[3] != "0"}
pop = {r[1]: r[5] for r in rows}
with open("phase_samples.txt", "w") as fh:
    for s in sorted(pop):
        if s not in hold:
            fh.write(s + "\n")
with open("rfmix_samples.tsv", "w") as fh:
    for g, pops in GROUPS.items():
        for p in pops:
            pick = sorted(s for s in pop if pop[s] == p and s not in hold and s not in children)[:PER_POP]
            assert len(pick) == PER_POP, (p, len(pick))
            for s in pick:
                fh.write(f"{s}\t{g}\n")
z = zipfile.ZipFile("maps.zip")
for c in range(1, 23):
    name = next(n for n in z.namelist() if n.endswith(f"plink.chr{c}.GRCh38.map"))
    lines = z.read(name).decode().splitlines()
    with open(f"plink.chr{c}.map", "w") as a, open(f"rfmix.chr{c}.map", "w") as b:
        for l in lines:
            ch, _, cm, pos = l.split()
            ch = ch if ch.startswith("chr") else f"chr{ch}"
            a.write(f"{ch}\t.\t{cm}\t{pos}\n")
            b.write(f"{ch}\t{pos}\t{cm}\n")
EOF
rm maps.zip
cut -f1 rfmix_samples.tsv > rfmix_ids.txt

one_chrom() {
    c=$1
    [ -s rfmix.chr$c.vcf.gz.tbi ] && return 0
    f=1kGP_high_coverage_Illumina.chr$c.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
    want=$(awk -v f="$f" '$1 == f { print $2 }' manifest.txt)
    for attempt in 1 2 3; do
        curl -sfL "$BASE/$f" | tee >(md5sum | cut -d' ' -f1 > chr$c.md5) \
            | bcftools view -m2 -M2 -v snps -i 'INFO/AF>=0.05 && INFO/AF<=0.95' -Ou \
            | bcftools annotate -x INFO,^FORMAT/GT -Ob -o common.chr$c.bcf || true
        sleep 2   # let the md5 subshell finish
        [ "$(cat chr$c.md5 2>/dev/null)" = "$want" ] && break
        echo "chr$c: download attempt $attempt failed (md5 $(cat chr$c.md5 2>/dev/null) != $want)" >&2
        [ "$attempt" = 3 ] && return 1
    done
    bcftools index -f common.chr$c.bcf
    bcftools query -f '%CHROM\t%POS\n' common.chr$c.bcf \
        | awk -v gap=3000 '$2 - last >= gap { print; last = $2 }' > sites.chr$c.txt
    bcftools view -T sites.chr$c.txt -S phase_samples.txt -Oz -o phase.chr$c.vcf.gz -W=tbi common.chr$c.bcf
    bcftools view -T sites.chr$c.txt -S rfmix_ids.txt -Oz -o rfmix.chr$c.vcf.gz -W=tbi common.chr$c.bcf
    rm -f common.chr$c.bcf* chr$c.md5
    echo "chr$c: $(wc -l < sites.chr$c.txt) sites" >&2
}
export -f one_chrom; export BASE
seq 1 22 | xargs -P "$THREADS" -I{} bash -c 'one_chrom {}'

# sites for reading the person's genotypes from the gVCF (pgs_engine.py target)
printf "#CHROM\tPOS\tID\tREF\tALT\n" > sites.pvar
for c in $(seq 1 22); do
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' rfmix.chr$c.vcf.gz \
        | awk 'BEGIN { OFS = "\t" } { sub(/^chr/, "", $1); print $1, $2, $1 ":" $2 ":" $3 ":" $4, $3, $4 }'
done >> sites.pvar
rm -f sites.chr*.txt
stamp "$(( $(wc -l < sites.pvar) - 1 )) sites, $(wc -l < phase_samples.txt) phasing genomes, $(wc -l < rfmix_samples.tsv) RFMix reference genomes"
