# Development and validation

How v2 was built and checked, written so it can be rerun. Paths are placeholders:

```bash
V2=/path/to/genomeportrait/v2          # this folder
SCRATCH=/fast/local/disk               # Nextflow work dirs and test outputs (not a FUSE/NFS mount)
REFBASE=$SCRATCH/refs                  # a reference base for tests
PANEL=/path/to/pgsc_HGDP+1kGP_v1       # extracted pgsc_calc reference panel (GRCh38 files)
```

Requirements on the test machine: Nextflow ≥ 24.04 with Java 17–21, Singularity/Apptainer,
`python3` with numpy, `rsync`. A plink2 Singularity image for the standalone engine
(`PLINK2_SIF`).

## Unit tests

```bash
python3 $V2/tests/test_pgs_target.py
```

The gVCF → reference-allele-space conversion on 17 hand-made cases (reference blocks,
zero depth, low GQ, uncalled, deletions spanning later records, indel spellings,
other-allele sites, missing chromosomes).

## Polygenic engine outside Nextflow, against pgsc_calc

```bash
python3 $V2/scripts/pgs_fetch.py --panel $V2/assets/pgs/panel_lifestyle.tsv --outdir $SCRATCH/scorefiles
```

```bash
PLINK2_SIF=/path/plink2.sif THREADS=16 bash $V2/tests/pgs_standalone.sh $SCRATCH/pgs $SCRATCH/scorefiles $PANEL <sample.g.vcf.gz> <sample>
```

```bash
python3 $V2/tests/compare_pgsc_calc.py $SCRATCH/pgs/<sample>.pgs_scores.tsv <pgsc_calc_results>/<sampleset>/score/<sampleset>_pgs.txt.gz
```

`pgs_standalone.sh` also accepts a pgsc_calc-style panel-space VCF (a file ending in
`target.vcf.gz`) in place of the gVCF, which isolates the engine from the conversion.
`TARGET_OPTS=--trust-all-blocks` reproduces how earlier targets treated zero-depth blocks,
for like-for-like comparisons. Results of the validation are in
[`polygenic.md`](polygenic.md).

Before a Nextflow run, compile the scripts with the report container's Python (3.11),
which rejects some f-strings the host's newer Python accepts:

```bash
singularity exec $REFBASE/containers/genomeportrait-report-1.0.1.sif sh -c 'cd /tmp/scripts && for f in *.py; do python3 -m py_compile $f || echo FAIL $f; done'
```

## Lint

```bash
cd $V2 && nextflow lint .
```

## Nextflow test runs

```bash
SCRATCH=$SCRATCH REFBASE=$REFBASE bash $V2/tests/run_test.sh smoke
```

FASTQ → bwa-mem2 → MarkDuplicates → CRAM → DeepVariant → genome facts → report, on the
nf-core chr21 test reads.

```bash
SCRATCH=$SCRATCH REFBASE=$REFBASE bash $V2/tests/run_test.sh modes
```

The smoke test's outputs fed back as two rows: alignment only (variant calling without
MarkDuplicates) and alignment + gVCF + VCF (interpretation only).

A portrait run on a real genome: write a samplesheet with `bam`, `gvcf` and `vcf`, and a
small config switching off what is not needed, then

```bash
SCRATCH=$SCRATCH REFBASE=$REFBASE bash $V2/tests/run_test.sh portrait --input samples.csv -c portrait.config --pgs_reference_panel $PANEL
```

## Seeding a reference base by hand

A `storeDir` is only honoured when every declared output is present, `versions.yml`
included. To reuse a FASTA or a container image that is already on disk:

```bash
mkdir -p $REFBASE/GRCh38/genome/genome $REFBASE/containers
```

```bash
cp Homo_sapiens_assembly38.fasta $REFBASE/GRCh38/genome/genome/ && printf '"DL_FASTA":\n    wget: seeded\n' > $REFBASE/GRCh38/genome/versions.yml
```

```bash
cp genomeportrait-report-1.0.1.sif $REFBASE/containers/ && touch $REFBASE/containers/report_container_1.0.1.ready && printf '"BUILD_REPORT_CONTAINER":\n    singularity: seeded\n' > $REFBASE/containers/versions.yml
```

## Personal data

Test runs on a real genome produce personal results. Keep them, and anything naming the
sample, out of this repository: the `local/` folder next to this file's parent is
git-ignored for exactly that.
