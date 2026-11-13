# genomeportrait: Changelog

## v2.0.0dev — 2026-11-12

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
