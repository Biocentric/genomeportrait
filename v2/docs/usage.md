# Usage

## Requirements

- Nextflow ≥ 24.04 with Java 17–21, Singularity/Apptainer (or Docker).
- ~60 GB for references with the defaults (GRCh38, VEP cache, ancestry panel, PGS reference
  panel), more with BQSR or SnpSift gnomAD annotation.
- Work directory on local disk. The pipeline itself can live anywhere, including a mount
  without exec bits: its scripts are staged into each task and run through `python3`.

## Samplesheet

```csv
sample,fastq_1,fastq_2,bam,gvcf,vcf,sex,platform
```

| Column | Meaning |
|---|---|
| `sample` | Identifier used in every filename, in read groups and in the report. Use a pseudonym if the outputs may be shared. |
| `fastq_1`, `fastq_2` | Gzipped reads. `fastq_2` empty for single-end or long reads. |
| `bam` | A `.bam` or `.cram`, coordinate-sorted, duplicate-marked and indexed (`.bai` / `.crai` beside it). |
| `gvcf` | DeepVariant or HaplotypeCaller gVCF (`.g.vcf.gz` + `.tbi`). Needs `bam` and `vcf`. |
| `vcf` | Normalised variant VCF (`.vcf.gz` + `.tbi`). Needs `bam` and `gvcf`. |
| `sex` | `XX`, `XY` or `unknown` (default). `unknown` is resolved from depth. |
| `platform` | `illumina` (default; also right for DNBSEQ/MGI), `ont`, `pacbio`. |

What runs depends on which columns are filled:

- **reads only** — everything.
- **`bam` only** — variant calling and interpretation. The alignment is used as it is.
  This differs from v1, which ran MarkDuplicates and BQSR again on a supplied BAM.
- **`bam` + `gvcf` + `vcf`** — interpretation only (portrait mode). Giving `gvcf` without
  `vcf`, or either without `bam`, is an error: polygenic scoring needs the gVCF's
  reference blocks, annotation needs the variant VCF, and STR, HLA, lineage and
  pharmacogenomics read the alignment.

To re-portray a genome that v2 already processed, point the three columns at
`<sample>/alignment/*.cram`, `<sample>/variants/*.g.vcf.gz` and `<sample>/variants/*.norm.vcf.gz`
in its results folder (v1 used `alignment/<sample>/`, `variants/deepvariant/<sample>/`).

## Running

```bash
nextflow run /path/to/genomeportrait/v2 -profile singularity \
    --input samplesheet.csv --outdir results --reference_base /data/genomeportrait_refs
```

`-profile production` is `singularity` plus a 19-core resource plan in which a failing
interpretation module cannot abort the run.

> **`--reference_base` can be set anywhere**: on the command line, in a `-params-file` or in
> a `-c` config. Container image paths are resolved from it when each task starts. The one
> exception is the Singularity image cache, which is fixed when `nextflow.config` is read:
> if you set `reference_base` in a `-c` config, also set `singularity.cacheDir` there, or
> export `NXF_SINGULARITY_CACHEDIR` (which always wins).

## Main parameters

### Choosing what runs

`--skip_qc`, `--skip_annotation`, `--skip_ancestry`, `--skip_haplogroups`, `--skip_str`,
`--skip_prs`, `--skip_hla`, `--skip_pharmcat`, `--skip_traits`, `--skip_telomere`,
`--skip_mito_heteroplasmy`, `--skip_multiqc`; `--call_sv false`, `--call_cnv false`.

Genome facts (depth, sex, mtDNA copies) always run: other modules depend on the result.

`--skip_blood_groups` turns off blood-group prediction. ABO, Kell, Duffy and Kidd come from
their defining variants, RhD from read depth over the RHD gene (the common RhD-negative
deletion); silencing variants in an intact RHD gene are not detected. Like CYP2D6, Kell needs
an ALT-aware alignment: KEL lies under the ALT contig chr7_KI270803v1_alt.

`--skip_kir` turns off KIR typing. T1K types the KIR genes from the reads of the KIR region
on chr19 and of every chr19 ALT contig: GRCh38's primary haplotype lacks several KIR genes
(KIR2DL5, KIR2DS3/5, KIR3DS1), so reads from them can only sit on ALT contigs. The report lists
which genes are present, their alleles, and whether the genome carries a B haplotype. It is
gene content from read evidence, not a clinical KIR typing; copy number is not reported.

`--skip_archaic` turns off the Neanderthal / Denisovan step. hmmix (Skov et al. 2018) looks
for stretches of the genome dense in variants that Africans do not carry, trains its model on
the person's own genome and labels each segment by the archaic genome it shares most variants
with. The genome is unphased, so a segment means at least one of the two copies is archaic
there; totals are compared with Skov's published calls for 1000 Genomes merged the same way
(`assets/archaic/1000g_reference.tsv`, built by `assets/archaic/make_reference.py`). The first
run downloads ~3.3 GB of helper files and archaic genomes into `<reference_base>/hmmix`
(`--hmmix_helper_record`, `--hmmix_archaic_record`: Zenodo record IDs). For an African genome
the outgroup removes most of its own variation, and the archaic share is small by design.

`--skip_painting` turns off chromosome painting. The person's genotypes at the reference SNPs
come from the gVCF (hom-ref sites count), Beagle phases them into two copies per chromosome
against the 1000 Genomes 30x panel, and RFMix assigns each stretch of each copy to one of four
reference groups (African, European, East Asian, South Asian; 20 genomes from each of five
unadmixed populations per group). Ancestry from outside those groups is painted as the nearest
of them, and the two copies are not maternal and paternal: each chromosome is phased on its own.
The first run streams the panel (~30 GB) into a ~1 GB reference in `<reference_base>/painting`;
a genome then takes ~25 min on 16 threads. The genomes in `assets/painting/holdout.tsv` are left
out of both references so that `tests/validate_painting.sh` can paint them.

`--skip_23andme` turns off the 23andMe-format export (`<sample>/export/`). Its site list is
built once from Illumina's GSA manifests, downloaded at run time (`--array_manifest_b37/b38`).

`--skip_cyp2d6` turns off CYP2D6 calling with Cyrius. Cyrius needs the reads at the
CYP2D6/CYP2D7 locus on the primary assembly. v2 maps ALT-aware, so its own alignments are
fine; an alignment supplied in the samplesheet that was made against GRCh38 with ALT contigs
but without ALT-aware mapping (all v1 alignments) has most of those reads on
`chr22_KI270879v1_alt`, and Cyrius then makes no call — PharmCAT reports CYP2D6 as not
called rather than guessing. Re-align such genomes from FASTQ to get a CYP2D6 call.

`--skip_contamination` turns off the VerifyBamID2 contamination estimate (about 20 min per
30x genome; `--contamination_panel` can point at the smaller 10k-marker panel shipped in the same
container).

### Alignment and calling

| Parameter | Default | |
|---|---|---|
| `--aligner` | `bwa-mem2` | `minimap2` for long reads |
| `--variant_caller` | `deepvariant` | `haplotypecaller`, `clair3` |
| `--markdup` | `samtools` | `samtools`: bwa-mem2 → fixmate → sort → markdup → CRAM in one stream, no intermediate BAMs. `gatk`: sorted BAM, then GATK MarkDuplicates (v1). Long reads always use the second path. |
| `--alt_aware` | `true` | Map ALT-aware with the reference's bwa `.alt` file (fetched once). Keeps reads at loci with ALT haplotypes — CYP2D6, HLA, KIR — on the primary assembly. |
| `--skip_bqsr` | `true` | `--skip_bqsr false` turns BQSR on. DeepVariant gains nothing from it. |
| `--skip_trimming` | `false` | |
| `--use_gpu` | `off` | `auto` / `on` for Parabricks; see the GPU notes in `nextflow.config` |

### Polygenic scores

| Parameter | Default | |
|---|---|---|
| `--pgs_panel` | `assets/pgs/panel_lifestyle.tsv` | Scores to compute and what each means |
| `--pgs_ids` | none | Extra PGS Catalog IDs, comma-separated. Reported with percentile and Z only. |
| `--reference_panel` | none | Directory with an extracted `pgsc_HGDP+1kGP_v1` panel, used for ancestry **and** polygenic scores. Without it the 16 GB archive is downloaded once. |
| `--pgs_min_overlap` | `0.75` | A score with fewer of its variants genotyped is not reported |
| `--pgs_min_gq` | `0` | Hom-ref calls below this genotype quality count as missing |
| `--pgs_keep_ambiguous` | `false` | Keep A/T and C/G SNPs |
| `--pgs_impute` | `auto` | Frequencies that fill missing genotypes: `auto` (most similar population when ≥ 90%), `population`, `pooled` |
| `--pgs_percentile` | `auto` | `population`, `adjusted`, or chosen by how cleanly the sample matches one population |

Method, validation and how to extend the panel: [`polygenic.md`](polygenic.md).

### Resources

`--max_cpus`, `--max_memory`, `--max_time` cap every task.

## What is cached

Everything under `--reference_base` is computed once and reused by every later run:
genome and indices, VEP cache, ancestry panel, container images, and under `pgs/` the
scoring files, the reference-panel subset, the score matching and the reference
distribution of every score. The `pgs/<panel>_<digest>/` directory is keyed by the list of
scores, so adding a score builds a new cache and leaves the old one valid.

## Tests

```bash
nf-test test                                                         # whole pipeline in stub mode, ~20 s
SCRATCH=/fast/dir bash tests/run_test.sh smoke                      # FASTQ path on tiny data
SCRATCH=/fast/dir bash tests/run_test.sh modes                      # aligned-only and portrait rows (after smoke)
bash tests/pgs_standalone.sh <workdir> <scorefiles> <panel> <gvcf> <sample>   # PGS engine outside Nextflow
python3 tests/compare_pgsc_calc.py <native.pgs_scores.tsv> <pgsc_calc_pgs.txt.gz>
```
