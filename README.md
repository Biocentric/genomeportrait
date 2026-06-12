# ![genomeportrait](docs/images/genomeportrait_logo_light.png)

[![GitHub Actions CI Status](https://github.com/your-org/genomeportrait/actions/workflows/ci.yml/badge.svg)](https://github.com/your-org/genomeportrait/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/your-org/genomeportrait/actions/workflows/linting.yml/badge.svg)](https://github.com/your-org/genomeportrait/actions/workflows/linting.yml)
[![nf-core template](https://img.shields.io/badge/nf--core%20template-3.x-23aa62.svg)](https://nf-co.re/)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**genomeportrait** is a bioinformatics pipeline for **personal whole-genome sequencing (WGS) analysis**. It takes raw FASTQ or aligned BAM/CRAM files from a single individual and produces a state-of-the-art, recreational-genomics "portrait" covering data quality, germline variants, ancestry and lineage, forensic STR profiles, polygenic trait scores, HLA immune type, pharmacogenomics, and a panel of well-known "fun" traits — all assembled into one interactive HTML report.

> [!IMPORTANT]
> **Scope & ethics.** genomeportrait is built for **recreational, ancestry, and educational** genomics. It deliberately **does not** generate clinical disease-risk predictions or diagnostic output. Variant annotations are reported for interpretation/context only. Do not use this pipeline for medical decision-making. See [`docs/ethics.md`](docs/ethics.md).

The pipeline is built using [Nextflow](https://www.nextflow.io), follows the [nf-core](https://nf-co.re) template and module conventions, and runs entirely inside **Singularity/Apptainer** containers for full reproducibility. All reference genomes, databases, and annotation resources are **downloaded automatically** on first run and cached for reuse. Every tool uses a public Bioconda/BioContainers image; the one bespoke report image is **built locally by the pipeline itself** (no container registry to host — see [docs/usage.md](docs/usage.md#containers-no-registry-hosting-required)).

### Quick test

```bash
# fast wiring smoke test (no large downloads)
nextflow run . -profile test,singularity --outdir results_smoke

# full flow on real 30x NA12878 reads, restricted to chr21
bash assets/prepare_test_chr21.sh
nextflow run . -profile test_30x,singularity --outdir results_test30x
```

## Pipeline summary

![genomeportrait metro map](docs/images/genomeportrait_metromap.png)

1. **Reference preparation** (`PREPARE_GENOME`) — download & index GRCh38, known-sites VCFs, VEP cache, gnomAD/ClinVar, 1000G+HGDP reference panel, PGS Catalog scoring files, STR catalogs, HLA & PharmCAT resources. Cached between runs.
2. **Read QC & trimming** — [`fastp`](https://github.com/OpenGene/fastp) + [`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
3. **Alignment** — [`bwa-mem2`](https://github.com/bwa-mem2/bwa-mem2) (short read) / [`minimap2`](https://github.com/lh3/minimap2) (long read), sort & index with [`samtools`](http://www.htslib.org/)
4. **Duplicate marking & recalibration** — [`GATK4 MarkDuplicates`](https://gatk.broadinstitute.org/) + `BaseRecalibrator`/`ApplyBQSR`
5. **Alignment QC** — [`mosdepth`](https://github.com/brentp/mosdepth), `samtools stats`, [`Qualimap`](http://qualimap.conesalab.org/)
6. **Small-variant calling** — [`DeepVariant`](https://github.com/google/deepvariant) (germline SNV/indel) → optional joint genotyping with [`GLnexus`](https://github.com/dnanexus-rnd/GLnexus)
7. **Structural & copy-number variants** — [`Manta`](https://github.com/Illumina/manta) + [`CNVkit`](https://cnvkit.readthedocs.io/)
8. **Variant annotation** — [`Ensembl VEP`](https://www.ensembl.org/vep) + [`SnpSift`](https://pcingola.github.io/SnpEff/) against gnomAD v4, ClinVar, dbNSFP
9. **Ancestry & admixture** — [`PLINK2`](https://www.cog-genomics.org/plink/2.0/) merge with 1000 Genomes + HGDP, PCA projection, [`ADMIXTURE`](https://dalexander.github.io/admixture/) supervised global ancestry, ROH/inbreeding & [`KING`](https://www.kingrelatedness.com/) kinship
10. **Maternal & paternal lineage** — mtDNA haplogroup with [`HaploGrep2`](https://haplogrep.i-med.ac.at/), Y-chromosome haplogroup with [`Yleaf`](https://github.com/genid/Yleaf)
11. **Forensic STR profile** — [`HipSTR`](https://hipstr-tool.github.io/HipSTR/) CODIS core loci + [`ExpansionHunter`](https://github.com/Illumina/ExpansionHunter) repeat expansions
12. **Polygenic scores** — PGS Catalog scoring files applied with a [`pgsc_calc`](https://pgsc-calc.readthedocs.io/)-style `PLINK2 --score` engine (educational traits incl. height, educational attainment/cognitive performance proxies)
13. **HLA typing** — [`T1K`](https://github.com/mourisl/T1K) (class I & II 4-digit)
14. **Pharmacogenomics** — [`PharmCAT`](https://pharmcat.org/) star-allele calling & drug-gene annotations
15. **"Popular" trait lookup** — curated, non-clinical SNP panel (eye/hair colour, lactase persistence, earwax/body-odour, bitter taste, caffeine & alcohol metabolism, muscle fibre type, photic sneeze, etc.)
16. **Extras** — telomere length ([`Telseq`](https://github.com/zd1/telseq)), mtDNA heteroplasmy ([`mutserve`](https://github.com/seppinho/mutserve))
17. **Reporting** — [`MultiQC`](https://multiqc.info/) + a bespoke single-file interactive **GenomePortrait** HTML report with all tables and plots.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet `samplesheet.csv`:

```csv
sample,fastq_1,fastq_2,bam,sex,platform
JOHN_DOE,/data/john_R1.fastq.gz,/data/john_R2.fastq.gz,,XY,illumina
```

…or start from an existing alignment:

```csv
sample,fastq_1,fastq_2,bam,sex,platform
JANE_DOE,,,/data/jane.cram,XX,illumina
```

Now run the pipeline:

```bash
nextflow run your-org/genomeportrait \
   -profile singularity \
   --input samplesheet.csv \
   --outdir results \
   --genome GRCh38 \
   --reference_base /data/genomeportrait_refs
```

> [!WARNING]
> Provide pipeline parameters via the CLI or a Nextflow `-params-file`. Custom config files (including those provided by the `-c` Nextflow option) can be used to provide any configuration _except for parameters_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits

genomeportrait was assembled following nf-core community conventions. It re-uses nf-core/modules wrappers where available and stands on the shoulders of the dozens of open-source tools listed in [`CITATIONS.md`](CITATIONS.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in [`CITATIONS.md`](CITATIONS.md). This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the MIT license. If you use it, please cite nf-core (Ewels _et al._, _Nat Biotechnol_ 2020) and the individual tools.
