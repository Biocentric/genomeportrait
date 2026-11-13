# Output

Everything is written under `--outdir`. References, indices and container images are not:
they live in `--reference_base`.

```
results/
├── report/<sample>.genomeportrait.html     the unified portrait, self-contained
├── polygenic/<sample>/
│   ├── <sample>.pgs_overview.html          the polygenic section as its own page
│   ├── <sample>.pgs_interpreted.tsv        every number on that page, one row per score
│   ├── <sample>.pgs_interpreted.json       the same, structured (consensus per trait)
│   ├── <sample>.pgs_scores.tsv             raw sums, Z and percentiles before interpretation
│   ├── <sample>.pgs_population.tsv         similarity to each reference population; projected PCs
│   ├── <sample>.pgs_imputation.tsv         which allele frequencies filled missing genotypes
│   └── <sample>.pgs_target_stats.tsv       how many scored positions the gVCF genotyped, and why not
├── qc/
│   ├── genome_facts/<sample>/              depth, breadth, sex chromosomes, mtDNA copies
│   ├── mosdepth/<sample>/  fastp/<sample>/  fastqc/<sample>/  qualimap/<sample>/
├── alignment/<sample>/                     CRAM + index, MarkDuplicates metrics
├── variants/deepvariant|manta|cnvpytor/<sample>/
├── annotation/<sample>/                    VEP-annotated VCF and table
├── ancestry/<sample>/                      PCA, admixture, ROH, kinship
├── lineage/<sample>/                       mtDNA and Y haplogroups
├── forensic_str/<sample>/                  CODIS profile, repeat expansions
├── hla/<sample>/  pharmacogenomics/<sample>/  traits/<sample>/  extras/
├── multiqc/multiqc_report.html
└── pipeline_info/                          timeline, trace, DAG, software versions
```

Not published, unlike v1: trimmed FASTQ, the duplicate-marked BAM, reference downloads and
indices, and intermediates such as the panel-space VCF used for scoring.

## `<sample>.pgs_interpreted.tsv`

| Column | Meaning |
|---|---|
| `domain`, `trait_concept` | Where the score sits in the report |
| `PGS`, `trait_reported`, `source` | PGS Catalog ID, the authors' trait label, first author and year |
| `higher_score_means` | What a higher aligned score predicts |
| `coding_flipped` | `yes` if the score's own coding ran the other way and was reversed |
| `direction_verified` | `yes` if the direction was confirmed by weight correlation against an anchor score |
| `aligned_percentile`, `aligned_Z` | After flipping. Z is ancestry-adjusted. |
| `verdict` | `favourable`, `unfavourable` or `neutral` (also when \|Z\| < 0.5) |
| `R2_eval` | Variance explained in the score's published validation, if usable |
| `expected_shift_SD`, `expected_shift_units` | √R² × Z, and the same in trait units where a population SD is known |
| `OR_or_HR_per_SD`, `relative_odds_vs_average` | Binary traits: published per-SD odds ratio, and OR^Z |
| `AUC`, `eval_ancestry`, `gwas_ancestry`, `european_trained` | Validation context |
| `n_variants`, `match_%` | Size of the score and the share of its variants genotyped in this sample |

## `<sample>.pgs_scores.tsv`

`status` is `ok`, `below_min_overlap`, `not_scored` (no variant matched the reference) or
`no_variance`. `SUM` is comparable to `pgsc_calc`'s `SUM`. `Z_norm1` adjusts the mean for
ancestry, `Z_norm2` the mean and the variance; `Z_MostSimilarPop` and
`percentile_MostSimilarPop` rank the sample within the closest reference population.
