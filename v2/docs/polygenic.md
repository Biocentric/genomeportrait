# Polygenic scores: method, validation, extending the panel

## What is reported

For each score in the panel:

| Quantity | Meaning |
|---|---|
| `SUM` | Σ weight × effect-allele dosage over the score's variants. Meaningless alone. |
| `Z_norm1` | `SUM` with its ancestry-dependent **mean** removed, in SDs of the reference residuals |
| `Z_norm2` | `SUM` with ancestry-dependent **mean and variance** removed. The Z shown in the report. |
| `Z_MostSimilarPop`, `percentile_MostSimilarPop` | Position among reference genomes of the most similar population |

These are the quantities `pgsc_calc` v2 reports, computed here by a native engine
(`scripts/pgs_engine.py` plus `plink2`) because `pgsc_calc` needed 6 h 45 min for one
genome on this panel, most of it in single-threaded reformatting and matching.

## Method

**Reference.** The `pgsc_calc` HGDP+1kGP panel (GRCh38), 3,942 genomes, 3,330 after
removing relatives with the KING cutoff list shipped with it.

**Sample-independent, cached per panel** under `reference_base/pgs/<panel>_<digest>/`:

1. `PGS_FETCH` — harmonised GRCh38 scoring files. Positions are read from `hm_chr` /
   `hm_pos`; `chr_name` / `chr_position` are the authors' original coordinates, often GRCh37.
2. `PGS_SITES` — union of scored autosomal positions (11.5 M for the default panel).
3. `PGS_REFERENCE` — the panel restricted to those positions; allele frequencies; PCA on
   common (MAF ≥ 5%), LD-pruned (r² < 0.1 in 1 Mb), biallelic SNPs outside three long-range
   LD regions (~96,000 SNPs). Reference genomes are then *projected* onto the PCs with the
   same `plink2 --score` command later used for the sample, so both sit on one scale.
4. `PGS_MATCH` — each scoring row is matched to a reference variant by position **and**
   alleles: effect/other = ALT/REF or REF/ALT, then the reverse-complement of either.
   Strand-ambiguous A/T and C/G SNPs are dropped. When an allele pair fits two reference
   variants (A/AAT is both the insertion A>AAT and the deletion AAT>A), the more common
   one is taken. Every weight is rewritten as a weight on the reference ALT allele; a row
   whose effect allele is the REF contributes `w × (2 − d_alt)`, so its weight becomes
   `−w` and `2w` goes into a per-score constant.
5. `PGS_SCORE_REFERENCE` — `plink2 --score` for every score over the reference genomes.

**Per sample:**

6. `PGS_TARGET` — one pass over the gVCF. For each reference variant:

   | gVCF at that position | Genotype emitted |
   |---|---|
   | reference block, `GT 0/0`, depth > 0, GQ ≥ `--pgs_min_gq` | `0/0` |
   | variant record carrying the reference ALT | its allele count |
   | variant record with other alleles only | `0/0` — zero copies of that ALT |
   | no record, zero depth, uncalled, or GQ below the threshold | missing |

   The variant-only VCF cannot be used: it omits hom-ref sites, and because effect alleles
   are often the reference allele the loss is directional. Expanding the gVCF with
   `bcftools convert --gvcf2vcf` is not enough either: it keeps the `<*>` allele, so no
   allele can match.
7. `PGS_SCORE` — `plink2 --score` with `--read-freq` from the reference, so a missing
   genotype contributes the reference mean dosage: exactly what that variant contributes on
   average to the distribution the sample is compared against.
8. `PGS_REPORT` —
   - *Most similar population*: Gaussian classifier on the first 5 PCs, one mean and
     covariance per reference super-population.
   - *`Z_norm1`*: ordinary least squares of `SUM` on the first 4 PCs in the reference;
     the sample's residual over the residual SD.
   - *`Z_norm2`*: a heteroscedastic normal model, mean linear in the PCs and log-variance
     linear in the PCs, fitted by maximum likelihood (Khan et al. 2022; the "full
     likelihood" model of `pgsc_calc`).
   - A score is dropped when fewer than `--pgs_min_overlap` (75%) of its variants are
     genotyped in the sample.

Autosomes only: chrX needs sex-aware ploidy and almost no score in the catalogue carries it.

## Interpretation

`scripts/pgs_interpret.py`, driven by two tables.

`assets/pgs/panel_lifestyle.tsv` — one row per score: `pgs_id`, `concept`, `reversed`,
`direction_verified`.

`assets/pgs/concepts.tsv` — one row per trait concept: domain, title, what higher and
lower mean, `good` (+1 higher is healthier, −1 lower is, 0 neither), `kind` (`q`
quantitative, `b` binary), `unit`, `population_sd`, `r2_max`, a note.

Rules (unchanged from `prs-selfscore`):

- **Direction.** PGS Catalog weights follow the source GWAS coding, which for UK Biobank
  ordinal fields can run against the label (alcohol frequency: 1 = daily … 6 = never).
  Seven scores are marked `reversed`; each was confirmed by correlating its SNP weights
  with an anchor score of unambiguous direction (r between −0.29 and −0.62).
- **Consensus.** Per concept, the median aligned Z over European-trained scores. Scores
  trained elsewhere are shown and left out.
- **Effect size, quantitative traits.** Expected shift = √R² × Z, R² from the score's
  published validation (European evaluations first). R² is not used when it comes from the
  score's own training cohort (`assets/pgs/r2_exclusions.tsv`), when it is reported as a
  percentage, or when it exceeds the concept's `r2_max` (a heritability ceiling).
- **Effect size, binary traits.** Relative odds = OR^Z, using only ORs between 1/3 and 3;
  larger ones are top-versus-bottom-decile contrasts, not per-SD.
- **Verdict.** Neutral when |Z| < 0.5 or the concept has no health direction.

The overview page adds notes it derives from the data: scores that were flipped, scores
that disagree with their trait's consensus by 1.5 SD or more and have no usable published validation (no R² or per-SD odds ratio), scores from
other ancestries, and scores that could not be reported.

## Validation against pgsc_calc

Reference run: `pgsc_calc` v2.3.0 on a ~21x DNBSEQ genome (2026-08), same 204 scores,
same reference panel; 196 scores passed its matching threshold. Reproduce with
`tests/pgs_standalone.sh` and `tests/compare_pgsc_calc.py`.

**A. Same genotypes.** The native engine run on the target VCF `pgsc_calc` itself used,
so any difference is the engine's (191 scores in both):

| Quantity | r | mean \|diff\| | max \|diff\| |
|---|---|---|---|
| `SUM` | 0.99999 | | |
| `Z_norm1` | 0.9996 | 0.023 | 0.154 |
| `Z_norm2` | 0.9996 | 0.021 | 0.147 |
| `Z_MostSimilarPop` | 0.9996 | 0.019 | 0.133 |
| `percentile_MostSimilarPop` | 0.9996 | 0.5 points | 4.0 points |

Matching statistics agree: 146.2 M of 156.3 M rows matched (`pgsc_calc` 145.9 M), 6.75 M
strand-ambiguous rows dropped (6.77 M), 39.6 thousand strand flips (39.0). Most similar
population EUR in both. Fed these Z values, the interpretation step reproduces
`prs-selfscore`'s interpreted table exactly in every categorical column (direction,
flipped, verdict, R², OR, ancestry flags) for all 191 scores.

**B. Full chain from the gVCF, like for like.** The whole v2 path (own gVCF conversion,
matching, scoring, adjustment), with one switch set to treat zero-depth blocks the way the
August target did (`--trust-all-blocks`, a comparison-only option):

| Quantity | r | mean \|diff\| | max \|diff\| |
|---|---|---|---|
| `SUM` | 0.99996 | | |
| `Z_norm1` | 0.9996 | 0.014 | 0.27 |
| `Z_norm2` | 0.9997 | 0.013 | 0.20 |
| `Z_MostSimilarPop` | 0.9997 | 0.011 | 0.24 |
| `percentile_MostSimilarPop` | 0.9998 | 0.3 points | 5.2 points |

All 196 `pgsc_calc` scores are reproduced. Five more that `pgsc_calc` dropped at 71–74%
overlap pass here (201 in total), because indels are recovered that the August conversion
skipped. The largest differences are in those small indel-rich scores.

**C. What v2 does differently by default, and what it changes.** DeepVariant writes
`GT 0/0` for reference blocks it has no reads for (`GQ 1`, `MIN_DP 0`). In this genome
2.0% of the scored positions (233,817 of 11.7 M) lie in such blocks. The August target
counted them as homozygous reference. That is not neutral: it replaces "unknown" with
"zero copies of the alternate allele" at 2% of every score, which moves a score by
whatever those sites are worth. v2 leaves them missing, so they take the reference mean.

Effect on `Z_norm2`, 197 scores: mean |shift| 0.14, 36 scores move by more than 0.25,
maximum 0.67; 17 scores cross the |Z| = 0.5 "around average" line. The shift is not
random: zero-depth sites sit disproportionately in GC-rich and repetitive regions, and
replacing "unknown" with "zero alternate alleles" pushes every score towards whatever the
reference allele is worth. At trait level, 13 of 41 consensus values moved by 0.25 SD or
more against the earlier run; none reversed direction. **Percentiles from a target built
the old way carry this bias.**

## Low-depth genomes and `--pgs_min_gq`

Zero-depth blocks are always left missing. `--pgs_min_gq` goes further and also drops
hom-ref calls whose genotype quality is below a threshold. Hom-ref evidence at the scored
positions of the 21x genome, by GQ of the reference block:

| GQ | sites | share |
|---|---|---|
| 0–9 | 26,697 | 0.3% |
| 10–19 | 40,018 | 0.5% |
| 20–39 | 478,289 | 6.0% |
| 40+ | 7,487,801 | 93.2% |

Raising the threshold changes little here: against the default, `--pgs_min_gq 10` moves
`Z_norm2` by 0.018 on average (maximum 0.15) and `--pgs_min_gq 20` by 0.025 (maximum 0.15);
no score moves by more than 0.25 or changes band. The default is therefore 0. For a genome
below ~10x, where thin blocks are a much larger share, 10 is a reasonable setting;
`<sample>.pgs_target_stats.tsv` reports the same GQ table for every run, so the share at
risk can be read before deciding.

A missed heterozygote is not biased in one direction by itself (it is called hom-ref or
hom-alt with similar probability), which is why thin evidence costs precision rather than
shifting scores, unlike the zero-depth case above.

## Extending the panel

**Add a score to an existing concept.** Append a row to the panel table:

```
PGS00XXXX	ldl	0	0	LDL cholesterol
```

**Add a new trait.** Add a row to `concepts.tsv`, then panel rows that refer to it. If the
trait is quantitative and you want effects in units, give `unit` and `population_sd`; give
`r2_max` a value near the trait's SNP heritability so an inflated R² is not used.

**Check the direction of a new score** when its trait is a questionnaire item or an
ordinal code. Correlate its weights with a score of known direction over shared variants;
a negative correlation means `reversed = 1`. `scripts/pgs_weight_corr.py` does
this.

**Try a score without curating it.** `--pgs_ids PGS00XXXX` scores it and reports
percentile and Z under "Other scores", with no direction or verdict.

Changing the list of scores changes the panel digest, so the sample-independent steps run
again for the new list (about 17 minutes); the previous cache stays valid for the old list.
Editing a concept, a note or a `reversed` flag does not invalidate anything.

## Limits

- **Missing genotypes take the pooled reference mean**, not the mean of the sample's own
  population. For this European genome, imputing the ~2.6% missing sites from European
  frequencies instead would shift `Z_norm2` by 0.04 on average (95th percentile 0.09,
  maximum 0.29). Population-specific imputation is on the v2.1 list.
- Percentiles are relative to the reference panel's closest super-population, which for
  most people of European ancestry means ~670 genomes from several European populations.
- A sample that is not a clean match to one super-population gets its percentiles from the
  ancestry-adjusted Z instead (`--pgs_percentile auto`); those are less certain.
- The projection uses plain `plink2 --score` on PC loadings, which shrinks a new sample's
  coordinates slightly compared with the reference genomes that defined the PCs.
  `pgsc_calc` corrects for this (FRAPOSA OADP); with ~96,000 SNPs and only the leading PCs
  in use, the agreement above shows the effect is negligible here.
- The population SDs used for unit conversion are approximate UK Biobank values.
- The cache for one panel takes about 13 GB under `reference_base/pgs/`.
- Validated on one genome of European ancestry. The engine has not been compared with
  `pgsc_calc` on an admixed or non-European genome.
