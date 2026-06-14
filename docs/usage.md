# your-org/genomeportrait: Usage

## Introduction

genomeportrait analyses a **single individual's** whole genome and produces a recreational/ancestry/educational "portrait". It accepts raw FASTQ or pre-aligned BAM/CRAM and runs entirely in Singularity containers.

## Samplesheet input

You must provide a CSV samplesheet with `--input`. One row per sample (you can run several individuals at once; each gets its own report).

```bash
--input '[path to samplesheet file]'
```

| Column     | Description |
|------------|-------------|
| `sample`   | Unique sample identifier (no spaces). |
| `fastq_1`  | Gzipped FASTQ, read 1 (paired-end) or the only file (long read). Optional if `bam` given. |
| `fastq_2`  | Gzipped FASTQ, read 2. Leave blank for single-end / long read. |
| `bam`      | Pre-aligned `.bam`/`.cram`. Provide this *instead of* FASTQs to skip alignment. |
| `sex`      | `XX`, `XY`, or `unknown`. Drives Y-haplogroup and STR ploidy handling. |
| `platform` | `illumina`, `ont`, or `pacbio`. Long-read platforms switch the pipeline to minimap2 + Clair3. |

Example:

```csv
sample,fastq_1,fastq_2,bam,sex,platform
JOHN_DOE,/data/john_R1.fastq.gz,/data/john_R2.fastq.gz,,XY,illumina
JANE_DOE,,,/data/jane.cram,XX,illumina
```

## References

On the first run, genomeportrait downloads every reference genome, annotation database,
and population panel into `--reference_base` (default `references/` next to the pipeline).
Subsequent runs re-use the cache (`storeDir`), so the multi-hundred-GB download only happens once.

To populate the cache once and then run fully offline:

```bash
# first run on a machine with internet — downloads everything
nextflow run your-org/genomeportrait -profile singularity --input s.csv --outdir out \
    --reference_base /shared/genomeportrait_refs

# later runs (same or air-gapped machine) re-use the cache
nextflow run your-org/genomeportrait -profile singularity --input s.csv --outdir out \
    --reference_base /shared/genomeportrait_refs --download_references false
```

> The ancestry/admixture module needs the large 1000G+HGDP panel and is **off by default**
> (`--skip_ancestry true`). Enable it with `--skip_ancestry false`.

## Running the pipeline

```bash
nextflow run your-org/genomeportrait \
    -profile singularity \
    --input samplesheet.csv \
    --outdir results \
    --genome GRCh38 \
    --skip_ancestry false
```

### Choosing which analyses run

Every interpretation layer has a `--skip_*` switch (see `nextflow_schema.json`), e.g.
`--skip_hla`, `--skip_str`, `--skip_prs`, `--skip_pharmcat`, `--skip_traits`,
`--skip_telomere`, `--skip_mito_heteroplasmy`.

### Long-read data

Set `platform` to `ont`/`pacbio` in the samplesheet and `--dv_model_type ONT_R104` (or `PACBIO`).
The pipeline aligns with minimap2 and calls variants with Clair3.

### Polygenic scores

`--pgs_ids` is a comma-separated list of [PGS Catalog](https://www.pgscatalog.org/) IDs.
Defaults point at educational/recreational traits (height, cognitive/educational-attainment
proxies, etc.). Scores are reported as a **relative ranking only**, never as absolute risk.

## Containers (no registry hosting required)

All third-party tools use public Bioconda/BioContainers images. The only bespoke image is the
small **report container** (pandas + matplotlib), which the pipeline **builds for you on the
first run** from [`containers/report/Singularity.def`](../containers/report/Singularity.def) and
caches under `--reference_base/containers/`. Nothing needs to be hosted anywhere.

- Under `-profile singularity`/`apptainer` it runs `apptainer build --fakeroot …` (falling back
  to `singularity build`). Unprivileged `--fakeroot` builds need a recent Apptainer; if your site
  forbids local builds, set `--build_report_container false` and either pre-build the `.sif`
  yourself (`apptainer build <path-from --report_sif> containers/report/Singularity.def`) or use
  `-profile conda`, which provisions the same environment with no image at all.
- Under `-profile docker` it runs `docker build` to create `genomeportrait-report:1.0.0` locally.

## Testing the pipeline

Two test profiles are provided:

```bash
# 1) Fast wiring smoke test (tiny chr21 reads; align → call → report; no big DB downloads)
nextflow run . -profile test,singularity --outdir results_smoke

# 2) Full-flow test on REAL 30x data restricted to chr21 (NA12878/HG001)
bash assets/prepare_test_chr21.sh          # one-off: streams + slices chr21, writes ./testdata
nextflow run . -profile test_30x,singularity --outdir results_test30x
```

`test_30x` runs essentially the whole pipeline (QC, alignment, DeepVariant, SV/CNV, VEP
annotation, haplogroups, STR, HLA, PharmCAT, traits, telomere, mtDNA, report). It uses the full
GRCh38 references, so the first run performs the one-time reference download (including the
~25 GB VEP cache). Genome-wide-only stages (ancestry, polygenic scores) are off because chr21-only
input isn't meaningful for them — flip `--skip_ancestry false --skip_prs false` on whole-genome data.

## GPU acceleration (NVIDIA Parabricks) — `dev` branch

When a Parabricks-capable GPU is available, the pipeline can swap the CPU alignment/calling
toolchain for [NVIDIA Parabricks](https://docs.nvidia.com/clara/parabricks/) — typically a
10–50× speed-up for `fq2bam` (bwa-mem + MarkDuplicates + BQSR) and `deepvariant`.

```bash
nextflow run . -profile gpu,singularity --input s.csv --outdir results        # auto-detect
nextflow run . -profile singularity --use_gpu on  ...                          # require a GPU
```

- `--use_gpu off` (default) — CPU tools only.
- `--use_gpu auto` — a `CHECK_GPU` preflight queries `nvidia-smi`; if **any** GPU meets
  `--gpu_min_memory_mb` (default 16000) **and** `--gpu_min_compute_cap` (default 7.0, i.e.
  Volta+), short-read samples run on Parabricks; otherwise it transparently falls back to CPU.
- `--use_gpu on` — same detection, but the run **fails fast** if no capable GPU is found.

Requirements: an NVIDIA GPU with **compute capability ≥ 7.0** and **≥ 16 GB VRAM** (T4, V100,
A100, L4, A6000, H100, …; note Pascal cards such as the **P40/P100 are *not* supported by
Parabricks 4.x**), the NVIDIA Container Toolkit, and Singularity/Apptainer started with `--nv`
(handled by the `gpu` profile). The Parabricks container is pulled from
`nvcr.io/nvidia/clara/clara-parabricks`. Long-read (ONT/PacBio) samples always use the CPU path.

## Core Nextflow arguments

- `-profile` — configuration profile. Use `singularity` (recommended) or `apptainer`, `docker`, `conda`.
- `-resume` — restart from cached results.
- `-c` — supply a custom config for compute/resources (not parameters).

## Reproducibility

Pin a release with `-r 1.0.0`. All tools are version-pinned containers, so a pinned
release + a populated reference cache is byte-for-byte reproducible.
