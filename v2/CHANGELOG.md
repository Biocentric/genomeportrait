# genomeportrait: Changelog

## v2.0.0dev — 2026-11-12

### Since the first v2 commit

- Fixed: boolean parameters given on the command line (`--skip_qc true`, `--skip_bqsr false`)
  failed validation ("Value is [string] but should be [boolean]"): Nextflow 25.10+ passes them as
  strings and nf-schema 2.2.0 did not convert them. nf-schema is now 2.8.0; the per-genome
  `references` table is excluded from the unrecognised-parameter check.
- Chromosome painting (`PAINT_CHROMOSOMES`, `PAINTING_REFERENCE`, `BUILD_PAINTING_CONTAINER`,
  `--skip_painting`): genotypes at ~763,000 common SNPs from the gVCF, phased with Beagle 5.4
  against the 1000 Genomes 30x panel, painted with RFMix v2 as African, European, East Asian or
  South Asian (400 unadmixed reference genomes). Report section with genome-wide shares, tract
  counts and lengths, and a karyogram of both copies of every chromosome. The reference is
  streamed once from the panel (~30 GB, MD5-checked) into ~1 GB.
- Neanderthal and Denisovan segments with hmmix (`ARCHAIC_HMMIX`, `HMMIX_RESOURCES`,
  `BUILD_HMMIX_CONTAINER`, `--skip_archaic`): the person's own variants, minus those seen in
  the African-rich outgroup, through hmmix create_ingroup, train and decode; segments labelled
  by the variants they share with the Altai, Vindija, Chagyrskaya and Denisova genomes. Report
  section with totals, the closest Neanderthal, the ten longest segments, a chromosome map and
  the 1000 Genomes populations with similar totals. Helper files and archaic genomes (~3.3 GB,
  Zenodo, CC-BY 4.0) are fetched once with MD5 checks.
- KIR gene content and alleles with T1K against IPD-KIR (`BUILD_T1K_KIR`, `BAM_KIR`,
  `--skip_kir`): reads from the KIR region on chr19 and all chr19 ALT contigs (several KIR
  genes, KIR2DL5 among them, are absent from the primary GRCh38 haplotype), a gene table with
  the A/A or B/x haplotype group in a report section. Calls with T1K quality 0 are dropped.
- Blood groups: ABO, Kell, Duffy, Kidd from their defining variants (read from the gVCF), RhD
  from RHD read depth (`RHD_DEPTH`, `BLOOD_GROUPS`, `--skip_blood_groups`), in a report section.
- ALT-aware alignment (`--alt_aware`, default on): bwa-mem2 gets the reference's `.alt` file, so
  reads at loci with ALT haplotypes stay on the primary assembly. v1 lacked it; on NA12878 chr22
  the CYP2D6 locus has 8,935 instead of 2,696 reads and Cyrius calls *3/*68+*4 instead of nothing.
- 23andMe-format export at the 639,615 SNPs of the Illumina GSA v3.0 (`ARRAY_PANEL`, `EXPORT_23ANDME`,
  `--skip_23andme`): GRCh37 positions, plus-strand bases, hom-ref sites from the gVCF, haploid
  chrX/Y/M for XY genomes.
- CYP2D6 star alleles from the alignment with Cyrius (image built locally; Cyrius is PyPI-only),
  passed to PharmCAT as an outside call; PharmCAT now runs its preprocessor and main program as two
  steps because `pharmcat_pipeline` takes no outside calls. Needs ALT-aware alignments (see TODO).
- Contamination check: VerifyBamID2 estimates the fraction of DNA from another person (FREEMIX),
  shown in "Genome at a glance" with a verdict (clean < 2%, check 2-5%, contaminated > 5%).
  `tests/validate_contamination.sh` mixes two samples at known fractions.
- Short reads go from bwa-mem2 through samtools fixmate, sort and markdup straight to a CRAM
  (`BWAMEM2_MARKDUP_CRAM`, `--markdup samtools`, the default). No sorted BAM or duplicate-marked
  BAM is written (~240 GB at 40x). `--markdup gatk` keeps the GATK MarkDuplicates path.
- Publishing is off by default and every result lands in one folder per sample:
  `<outdir>/<sample>/{alignment,variants,qc,annotation,ancestry,lineage,str,polygenic,hla,
  pharmacogenomics,traits,extras}/` plus `<outdir>/<sample>/<sample>.genomeportrait.html`;
  `multiqc/` and `pipeline_info/` stay at the top. `tests/check_results_tree.py` checks it.
- Ancestry inputs fixed (both v1 bugs): hom-ref sites were called with ALT `.` and dropped at the
  merge (only 19,663 of 50,000 SNPs used, mostly variant sites), and the 50,000 SNPs were the first
  lines of the pruned list (chromosomes 1-3 only). Sites are now spread genome-wide and genotyped for
  the panel's alleles; the sample is named after the samplesheet, not the CRAM read group.
  Test genome: 99.2% EUR instead of 80%.
- Ancestry uses the same HGDP+1kGP reference panel as the polygenic scores (`--reference_panel`,
  prepared once in `PREPARE_GENOME`), with relatives removed; the separate 1000G+HGDP download
  (`DOWNLOAD_KGP_HGDP`) is gone. `--pgs_reference_panel(_url)` became `--reference_panel(_url)`.
  Hold-out check: `tests/validate_ancestry_holdout.sh`.
- Missing genotypes in polygenic scoring are filled from the most similar reference
  population's allele frequencies (`PGS_POP_FREQ`, `--pgs_impute`) instead of the pooled panel's.
  Removes a systematic error of up to 0.35 SD; `tests/validate_pgs_imputation.py` measures it.
- The sample is scored in one pass in Python (`pgs_engine.py score`) instead of one plink2 run per
  score: 40 s instead of 4 min 50 s for 204 scores. `PGS_MATCH` also saves each score as row
  numbers and weights (`<pgs>.npz`, cache directory `match_v2`); `PGS_SCORE` became
  `PGS_PROJECT` (PCA projection only). Percentile ties are compared with a tolerance.
- gVCF conversion runs one chromosome per CPU, reading the gVCF's tabix index directly (no
  external tool): 1 min 18 s instead of 5 min 49 s on a 21x genome, identical output.
- `report_sif`, `hipstr_sif` and `yleaf_sif` default to null and are resolved from `reference_base`
  when each task starts, so `reference_base` may be set in a `-c` config. The Singularity cache
  honours `NXF_SINGULARITY_CACHEDIR`.
- Report notes: a score that disagrees with its trait is described as having "no usable published
  validation (no R² or per-SD odds ratio)", and as evaluated only in non-European samples where
  that applies, instead of "no published validation" — most such scores do have performance records.

Next version of the pipeline, forked from v1 branch `feat/parabricks-pascal` @ `f4eb0ae`
and kept in its own folder (`v2/`). See [`PLAN.md`](PLAN.md) for the reasoning and
[`TODO.md`](TODO.md) for what is open.

### Added

- **Portrait mode.** `gvcf` and `vcf` samplesheet columns. A row with `bam` + `gvcf` + `vcf`
  runs interpretation only.
- **Polygenic scores, rebuilt** (`VCF_POLYGENIC`, `scripts/pgs_*.py`, `assets/pgs/`):
  - 204-score `lifestyle` panel as an explicit, versioned table, with a trait concept per score.
  - Native engine: allele-aware matching to the HGDP+1kGP reference, the sample rebuilt in
    the reference's allele space from the gVCF, `plink2 --score`, PCA projection,
    ancestry-adjusted `Z_norm1` / `Z_norm2`, percentile in the most similar population.
    Agrees with `pgsc_calc` v2.3.0 at r = 0.9997 over the full chain (see `docs/polygenic.md`).
  - Positions the gVCF has no reads for are left missing and take the reference mean. The
    2026-08 `pgsc_calc` target counted them as hom-ref; correcting that moves 36 of 197
    scores by more than 0.25 SD on the test genome.
  - Interpretation: flipped codings, per-trait consensus, expected effect size in units,
    verdict; from `prs-selfscore/interpret/`, now table-driven.
  - Overview page (`<sample>.pgs_overview.html`) and the same content inside the unified report.
  - Sample-independent work cached per panel under `reference_base/pgs/`.
  - PGS Catalog metadata snapshot shipped in `assets/pgs/catalog/`.
- **Genome facts** (`GENOME_FACTS`): depth, breadth, inferred sex chromosomes, mtDNA copies
  per cell. The resolved sex drives ExpansionHunter ploidy and whether Yleaf runs.
- `SAMTOOLS_BAM2CRAM`: the alignment deliverable is a CRAM whether or not BQSR runs.
- `tests/`: standalone PGS run, comparison against `pgsc_calc`, gVCF-conversion unit test, generic test launcher.

### Changed

- **A samplesheet `bam` is taken as analysis-ready.** It is called as it is; v1 ran
  MarkDuplicates and BQSR on it again.
- **BQSR is off by default** (`--skip_bqsr false` enables it). Known-sites files are only
  downloaded when something reads them.
- **One sort.** `SAMTOOLS_SORT` after the aligner is gone; both aligners already sort.
- The bwa-mem2 index is only built when a samplesheet row starts from FASTQ.
- `BAM_EXTRAS` converts CRAM to BAM only when telomere or heteroplasmy analysis is on.
- **Publishing**: references, indices, containers, trimmed FASTQ, the duplicate-marked BAM
  and step-to-step intermediates are no longer copied into `results/`.
- **Scripts are staged process inputs** in `scripts/`, run as `python3 <script>`. `bin/` is gone.
- Unified report restyled (light and dark, card layout), sections reordered.
- mosdepth moved out of `BAM_QC`; it always runs.
- Ethics statement updated: the default polygenic panel now includes health-related traits.
- Named `Biocentric/genomeportrait`; template placeholders removed.

### Fixed

- The report's variant overview section never rendered in v1: it looked for `variant_stats`
  in the filename while the file is `<sample>.bcftools_stats.txt`.
- PGS genotyping kept the gVCF's `<*>` allele, so alternate-allele effect weights at
  hom-ref sites were dropped rather than counted as zero copies.
- ExpansionHunter assumed `female` for an undeclared sex; Yleaf ran on every sample whose
  sex was undeclared.
- The index of a samplesheet CRAM was looked up as `<file>.bai`.

### Removed

- `PGSCATALOG_DOWNLOAD`, `PGS_GENOTYPE`, `PLINK2_VCF2PGEN`, `PLINK2_SCORE`,
  `PLINK2_SCORE_REFERENCE`, `PRS_REPORT`, `SAMTOOLS_SORT`, `bin/`.
- `--pgs_ids` no longer defines the panel; it adds scores to it.

## v1

See the changelog one level up.
