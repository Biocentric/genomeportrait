# Ethics, scope and responsible use

genomeportrait is a recreational, ancestry and educational genomics pipeline. Read this
before running it on anyone's data.

## What it is for

- Exploring your own genome out of curiosity: ancestry, lineage, traits, HLA type, how
  your inherited tendencies compare with people of similar ancestry.
- Learning how whole-genome analysis works from reads to report.

## What it is not

- **Not a medical or diagnostic device.** No part of the output may be used for diagnosis,
  screening, treatment or any clinical decision.
- The pharmacogenomics (PharmCAT) and annotation (ClinVar) sections are context. Any
  medication decision needs a clinician and an accredited clinical test.

## What changed in v2: health-related polygenic scores

v1 stated that it produced no disease-risk scores. **That is no longer true of the default
panel.** The v2 `lifestyle` panel includes scores for LDL, HDL and triglycerides, BMI,
resting heart rate, coronary artery disease, type 2 diabetes, metabolic syndrome, sleep and
frailty, and the report says which direction is healthier and how large the expected
effect is.

This was a deliberate choice, taken from the `prs-selfscore` work this panel comes from.
What keeps it within "educational" rather than "clinical":

- Every result is an **inherited tendency relative to people of similar ancestry**, stated
  as a shift from average, never as an absolute risk or probability.
- An effect size is only shown where the PGS Catalog holds a usable published validation;
  otherwise the report says "direction only". For most traits the score explains under
  10% of the variation between people, and the report says so on each card.
- A measured value (an actual LDL, HbA1c, heart rate) always outranks a score.
- Scores transfer poorly across ancestries. Scores trained on other ancestries than the
  sample's closest reference population are shown but kept out of every consensus.

If you want the v1 position back, point `--pgs_panel` at a panel without those rows; the
panel is a plain TSV. Monogenic findings (ACMG secondary findings, APOE, carrier status)
remain out of scope.

## Consent and privacy

- Only analyse genomes you own or have explicit, informed consent for.
- A genome identifies a person and reveals information about relatives. The forensic STR
  (CODIS) and kinship modules are identity-grade.
- Outputs are sensitive personal data. Keep `--outdir` and `--reference_base` private and
  do not upload reports to third-party services without consent.
- Inferred sex chromosomes are reported only as XX, XY or undetermined. The pipeline does
  not name sex-chromosome aneuploidies: that is a clinical statement, and a depth ratio
  from one short-read run does not support it.

## Trait panel

`assets/trait_snps.tsv` is limited to non-clinical characteristics (appearance, taste,
simple physiology). Behavioural associations in it are small, debated, and conversation
pieces only.
