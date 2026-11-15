# genomeportrait v2

A Nextflow pipeline that turns one person's whole-genome sequencing data into a single
HTML portrait: ancestry, maternal and paternal lineage, forensic STR profile,
ancestry-adjusted polygenic scores, HLA type, pharmacogenomics and a panel of well-known
traits.

This folder is the next version of the pipeline. It is self-contained and does not touch
the v1 files one level up. Why it exists and what is still open:
[`PLAN.md`](PLAN.md), [`TODO.md`](TODO.md), [`CHANGELOG.md`](CHANGELOG.md).

> **Scope.** Recreational, ancestry and educational genomics. Nothing in the report is a
> diagnosis. Polygenic scores are inherited tendencies relative to people of similar
> ancestry, not measurements. See [`docs/ethics.md`](docs/ethics.md).

## What is new in v2

| | v1 | v2 |
|---|---|---|
| Re-running interpretation on a finished genome | from BAM, repeating MarkDuplicates, BQSR and DeepVariant | **portrait mode**: give CRAM + gVCF + VCF, interpretation only |
| Polygenic scores | 4 scores, pooled percentile | **204-score panel**, ancestry-adjusted Z, direction, consensus per trait, effect size in units |
| PGS run time for one genome | n/a (`pgsc_calc` on the same panel: 6 h 45 min) | ~20 min for the first genome on a panel, ~3 min after that |
| Sample sex | samplesheet only; `unknown` treated as female for STRs | **inferred from depth**, used when the samplesheet does not say |
| Alignment | sorted twice, BQSR always, two whole-genome BAMs on disk | aligned, sorted and duplicate-marked straight to CRAM in one stream; BQSR opt-in |
| `results/` | also holds references, trimmed FASTQ, markdup BAM (~360 GB extra at 40x) | results only |
| Running from a mount without exec bits | no (`bin/` on `PATH`) | yes (scripts are staged inputs) |
| Report | dark table dump; variant overview section never rendered | light/dark, card layout, polygenic overview page |

## Quick start

```bash
nextflow run /path/to/genomeportrait/v2 -profile singularity \
    --input samplesheet.csv --outdir results --reference_base /data/genomeportrait_refs
```

The samplesheet decides, per row, how much of the pipeline runs:

```csv
sample,fastq_1,fastq_2,bam,gvcf,vcf,sex,platform
A,/data/a_R1.fastq.gz,/data/a_R2.fastq.gz,,,,unknown,illumina
B,,,/data/b.cram,,,XX,illumina
C,,,/data/c.cram,/data/c.g.vcf.gz,/data/c.vcf.gz,unknown,illumina
```

| Row | Given | Runs |
|---|---|---|
| A | reads | trim, align, mark duplicates, call variants, interpret |
| B | alignment | call variants, interpret (the alignment is taken as sorted and duplicate-marked) |
| C | alignment + gVCF + VCF | interpret only — **portrait mode** |

`sex` may be `unknown`; it is then inferred from chrX and chrY depth.

Full parameter and output descriptions: [`docs/usage.md`](docs/usage.md),
[`docs/output.md`](docs/output.md). How the polygenic scores are computed and validated:
[`docs/polygenic.md`](docs/polygenic.md).

## Pipeline

1. **Reference preparation** — downloads and indexes what the samplesheet needs, once,
   into `--reference_base`. No bwa-mem2 index is built for runs that start from alignments;
   no BQSR known-sites are fetched unless BQSR is switched on.
2. **Reads to calls** (rows with FASTQ) — fastp, FastQC, bwa-mem2 streamed through samtools
   fixmate/sort/markdup to CRAM (or minimap2 + GATK MarkDuplicates), DeepVariant (or HaplotypeCaller, Clair3). Parabricks on a capable GPU.
3. **Genome facts** — mosdepth; depth, breadth, sex chromosomes, mtDNA copies per cell.
4. **Interpretation** — ancestry (PCA, ADMIXTURE, ROH, KING), mtDNA and Y haplogroups,
   CODIS STR profile and repeat expansions, **polygenic scores**, HLA (T1K),
   pharmacogenomics (PharmCAT, with CYP2D6 from Cyrius), trait SNP panel, telomere length, mtDNA heteroplasmy,
   SV/CNV, VEP annotation.
5. **Report** — MultiQC plus one self-contained `<sample>.genomeportrait.html`, and the
   polygenic section on its own as `<sample>.pgs_overview.html`.

## Polygenic scores in one paragraph

The sample is rebuilt from its gVCF in the allele space of the HGDP+1kGP reference panel
(reference blocks become hom-ref, uncovered sites are left missing and take the reference
mean). Each PGS Catalog score is matched to the panel by position and alleles and applied
with `plink2 --score` to the sample and to 3,330 unrelated reference genomes. The score is
then adjusted for ancestry by modelling its mean and variance along the reference's
principal components. On the genome used for testing, the result agrees with `pgsc_calc`
v2.3.0 at r = 0.9997 across 196 scores when both treat the data the same way. By default
v2 treats one thing differently: positions the gVCF has no reads for are left missing
instead of being counted as hom-ref, which moves a number of scores by 0.25–0.7 SD
(details in [`docs/polygenic.md`](docs/polygenic.md)). Each score is finally given a
direction, a verdict and, where the PGS Catalog holds a usable validation, an expected
effect size.

## Layout

```
main.nf, workflows/, subworkflows/, modules/   the pipeline
scripts/                                       Python used by the pipeline (staged, never on PATH)
assets/pgs/                                    score panel, trait concepts, PGS Catalog snapshot
conf/                                          base, modules (publishing), production, test profiles
tests/                                         standalone PGS run, portrait-mode test, validation
docs/                                          usage, output, polygenic method, ethics, command log
```

## Status

v2.0.0dev. Tested here:

- the polygenic engine against `pgsc_calc`, same genotypes and full chain;
- the gVCF conversion on hand-made cases (`tests/test_pgs_target.py`);
- a portrait-mode run on a real 21x genome (genome facts, polygenic scores, report);
- the `test` profile from FASTQ on tiny data (align, mark duplicates, CRAM, DeepVariant, report);
- all three samplesheet modes in one run on tiny data.

Not yet done: a whole-genome run from FASTQ, and a portrait run with every interpretation
module switched on. The alignment-reading modules (ancestry, lineage, STR, HLA,
pharmacogenomics, traits, SV/CNV, annotation) are v1 code and were not exercised in v2's
tests. See [`TODO.md`](TODO.md).
