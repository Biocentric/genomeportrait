# genomeportrait v2 — plan

Written 2026-11-12. This folder (`genomeportrait/v2/`) is the next version of the
pipeline. It is self-contained: nothing in it is imported from, or writes to, the v1
files one level up. v1 keeps working exactly as it does today.

Progress is tracked in [`TODO.md`](TODO.md). Every command that was run to build and
test v2 is logged in [`docs/development.md`](docs/development.md).

---

## 1. Where v1 stands

Three copies of v1 exist and they are not the same code:

| Copy | State |
|---|---|
| `genomeportrait/` (this Filen folder) | `main` @ `733987b` plus ~45 uncommitted, CRLF-converted files. Predates the Parabricks work. |
| GitHub `Biocentric/genomeportrait` `dev` | `34e4dbd` — what this folder last fetched |
| NVMe clone used for the 42x run | branch `feat/parabricks-pascal` @ `f4eb0ae`, 11 commits ahead of `dev` — the only copy with the input-truncation fixes |

**v2 is forked from `f4eb0ae`**, the newest and the only one that has survived a
full 42x genome.

What works well in v1 and is kept: the reference auto-download with `storeDir`
caching, the GPU/CPU routing, the STR recovery logic, PharmCAT genotyping from the
alignment, and the per-module tolerance in the production profile.

What costs time or disk, or is simply wrong:

| # | Problem | Cost |
|---|---|---|
| 1 | Any change to an interpretation step means re-running from FASTQ or BAM: the samplesheet accepts a BAM but still re-runs MarkDuplicates, BQSR and DeepVariant on it | days per iteration at 42x |
| 2 | `BWAMEM2_MEM(sort=true)` followed by `SAMTOOLS_SORT` | one redundant whole-genome sort, ~120 GB |
| 3 | BQSR runs by default although DeepVariant neither needs nor wants it | one extra whole-genome pass, ~120 GB |
| 4 | Catch-all `publishDir` copies references, trimmed FASTQ and the markdup BAM into `results/` | ~360 GB per genome |
| 5 | Polygenic scores: 4 scores, percentile against a pooled panel with no ancestry adjustment, no direction, no effect size | the section says little |
| 6 | PGS genotyping keeps the gVCF `<*>` allele, so ALT-effect alleles at hom-ref sites are dropped instead of counted as 0 | biased coverage figure |
| 7 | Sample sex comes only from the samplesheet; `unknown` silently makes ExpansionHunter assume female and always runs Yleaf | wrong ploidy on chrX STRs |
| 8 | Report scripts are found through `bin/` on `PATH`, which needs the exec bit | cannot run from the Filen mount |
| 9 | Template leftovers (`your-org`, placeholder badges, a metro map that does not exist) | noise |

---

## 2. What v2 changes

### 2.1 Structure

1. **Portrait mode.** The samplesheet gains `gvcf` and `vcf` columns next to `bam`. A row
   that supplies all three skips trimming, alignment, duplicate marking and variant
   calling and goes straight to interpretation. A finished genome can be re-portrayed
   in about an hour instead of days, which is what makes iterating on the report
   practical. A row with `bam` alone is called as it is; v1 marked duplicates again.
2. **One sort, no BQSR by default.** `SAMTOOLS_SORT` after bwa-mem2 is removed.
   BQSR is off unless `--skip_bqsr false`; asking for it with DeepVariant is honoured
   with a warning. Without BQSR the duplicate-marked BAM is converted to CRAM as the
   deliverable. Known-sites files and the bwa-mem2 index are only prepared when
   something will use them.
3. **Publishing excludes what is not a result.** References, containers, indices,
   trimmed FASTQ, the duplicate-marked BAM and step-to-step intermediates are never
   copied to `results/`. (The default per-process rule is kept for everything else;
   a strict opt-in layout is listed under v2.1.)
4. **Scripts are staged as inputs** and run as `python3 script.py`. No exec bit, no
   `bin/` on `PATH`; the pipeline runs from a FUSE mount or a fresh clone alike.
5. **De-templated.** Named `Biocentric/genomeportrait`, version `2.0.0dev`.

### 2.2 Polygenic scores (the main upgrade)

Taken from `prs-selfscore/`: the score panel, the interpretation rules and the page
layout. Not taken: running `pgsc_calc` itself, which needed 6 h 45 min for one
sample, most of it in single-threaded Python that reformats scoring files.

v2 computes the same quantities natively:

| Step | Process | Cached? |
|---|---|---|
| Download harmonised GRCh38 scoring files | `PGS_FETCH` | yes, per panel |
| Union of scored positions | `PGS_SITES` | yes, per panel |
| Subset the HGDP+1kGP reference to those positions, drop relatives, PCA | `PGS_REFERENCE` | yes, per panel |
| Match every score to the reference by `chr:pos:ref:alt` (ambiguous A/T, C/G SNPs dropped, strand flips resolved) | `PGS_MATCH` | yes, per panel |
| Score the 3,330 reference genomes | `PGS_SCORE_REFERENCE` | yes, per panel |
| Rebuild the sample in the reference's allele space from the gVCF (hom-ref blocks become `0/0`, uncovered sites stay missing and are mean-imputed from reference frequencies) | `PGS_TARGET` | per sample |
| Score the sample, project it onto the reference PCs | `PGS_SCORE` | per sample |
| Ancestry-adjusted Z, percentile in the most similar population; direction, consensus per trait, effect size, verdict; overview page | `PGS_REPORT` | per sample |

Everything sample-independent is computed once and kept under
`reference_base/pgs/<panel>/`. The first genome pays for it; later genomes take
minutes.

**Panel.** `assets/pgs/panel_lifestyle.tsv` lists the 204 scores of the
`prs-selfscore` run explicitly (the run resolved them from 20 EFO terms plus 11
hand-picked IDs; an explicit list is reproducible, an ontology query is not). Each
row carries its trait concept and whether its coding is reversed.
`assets/pgs/concepts.tsv` holds the 42 trait concepts: domain, what "higher" means,
which direction is healthier, unit, population SD, heritability ceiling for R².
These were hard-coded in `interpret.py`; as tables they can be extended without
touching code. `--pgs_ids` adds extra scores, which are reported with percentiles
only until they are given a concept.

**Interpretation rules, unchanged from `prs-selfscore`:**

- Scores with reversed coding (7, verified by weight correlation against an anchor)
  are flipped so every bar reads the same way.
- Consensus per trait is the median aligned Z over European-trained scores.
- Continuous traits: expected shift = √R² × Z, converted to units where a population
  SD is known. R² is rejected when it is in-sample (Jung 2024) or above the trait's
  heritability ceiling.
- Binary traits: relative odds = OR<sup>Z</sup>, ORs outside 1/3–3 discarded.
- Catalog metadata and performance records ship as a snapshot in
  `assets/pgs/catalog/`, so interpretation works offline and is reproducible.

**Validation.** The earlier `pgsc_calc` run left its target VCF, reference panel and
all 196 results. The native engine is run on the same genome and
compared score by score. Acceptance: Pearson r ≥ 0.98 between native and `pgsc_calc`
Z across scores, and no trait changes verdict.

### 2.3 New analyses

Chosen for being cheap, testable on data already on disk, and useful to every other
module:

- **Genome facts** (`GENOME_FACTS`): sex karyotype inferred from X/Y/autosome depth,
  mtDNA copy number, mean depth and breadth. The inferred sex replaces `unknown` for
  ExpansionHunter and decides whether Yleaf runs.

Investigated and deferred (see §4): CYP2D6 calling, blood groups, archaic
introgression, chromosome painting, array-format export, contamination check.

### 2.4 Report

One self-contained HTML in the `prs-selfscore/overview.html` design: light and dark,
readable column width, cards. The polygenic section is the overview page itself
(diverging bars, verdict chips, per-score tables). Other sections keep their v1
content under the new styling. The polygenic section is also written standalone as
`<sample>.pgs_overview.html`, with every number in `<sample>.pgs_interpreted.tsv`.

---

## 3. Order of work

1. Fork `f4eb0ae` into `v2/`, apply §2.1.
2. PGS engine as standalone scripts; validate against `pgsc_calc` (§2.2).
3. Wrap the engine as Nextflow modules; port interpretation and overview page.
4. `GENOME_FACTS`.
5. Report restyle.
6. Nextflow tests: portrait mode on the real gVCF/CRAM; `-profile test` smoke test.
7. README, usage, command log.

---

## 4. Candidate analyses considered

| Analysis | Verdict | Why |
|---|---|---|
| Ancestry-adjusted PGS panel with interpretation | **v2.0** | The stated upgrade |
| Sex inference, mtDNA copy number | **v2.0** | Free from mosdepth; fixes a real v1 defect |
| CYP2D6 star alleles (Cyrius) fed to PharmCAT | v2.1 | PharmCAT cannot call CYP2D6 from a VCF; it is the most-used pharmacogene. Needs a container check and a truth sample. |
| Blood groups (ABO, RhD, Kell, Duffy, Kidd) | v2.1 | Popular and cheap, but ABO is minus-strand and RhD is a deletion call; needs validating against reference-panel phenotype frequencies before it is shown to anyone |
| One reference panel for ancestry and PGS | v2.1 | The PGS panel (16 GB) could replace the separate ancestry panel download. Changes validated v1 ancestry results, so it is its own step. |
| Streamed align → sort → markdup → CRAM | v2.1 | Removes two more ~120 GB intermediates; cannot be verified without a full-genome run |
| Contamination estimate (VerifyBamID2) | v2.1 | Cheap QC that guards every downstream call |
| 23andMe-format export at array positions | v2.2 | Lets the genome go into third-party tools; needs an array manifest and GRCh37 liftover |
| Neanderthal / Denisovan share | v2.2 | Popular; needs archaic reference genomes and a validated method |
| Chromosome painting (local ancestry) | v2.2 | Needs phasing against a haplotype panel; heavy |
| KIR typing (T1K) | v2.2 | T1K already in the pipeline; low demand |
| ACMG secondary findings, APOE, carrier status | not planned | Clinical. v1's ethics statement excludes it; a decision for you, not a default |

---

## 5. Risks, and how they turned out

- **The native PGS engine is new code.** Compared score by score with `pgsc_calc` on the
  same genome: r = 0.9997 over the full chain. The comparison also found a bug in the new
  conversion (fixed, unit-tested) and a bias in the August target (zero-depth blocks
  counted as hom-ref), which v2 does not reproduce. See `docs/polygenic.md`.
- **Low-depth genomes.** Measured on the 21x genome: the GQ threshold moves no score by
  more than 0.15 SD, so the default leaves it off. Not measured on a genome below 10x.
- **No full-genome run of v2 from FASTQ.** The read-to-call half ran on the tiny `test`
  profile only. Still open.
- **Interpretation modules other than polygenic scores were not run in v2.** They are v1
  code with a changed call signature; lint passes, but only a full portrait run proves it.
  Still open: the reference set they need sits on the RAID, which was not mounted.
- **Validated on one European genome.** Ancestry adjustment for a non-European or admixed
  genome is untested.
