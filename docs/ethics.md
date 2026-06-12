# Ethics, scope & responsible use

genomeportrait is a **recreational, ancestry, and educational** genomics pipeline. Please
read this before running it on anyone's data.

## What this pipeline is for

- Curiosity-driven exploration of your own genome (eye colour, lactase persistence, bitter
  taste, ancestry composition, haplogroups, HLA type, etc.).
- Learning how modern WGS analysis works end-to-end.
- Ancestry and deep-lineage (mtDNA / Y) inference.

## What this pipeline is **not** for

- **Not a medical or diagnostic device.** It does not produce disease-risk scores and must
  not be used for diagnosis, screening, treatment, or any clinical decision.
- The pharmacogenomics (PharmCAT) and variant-annotation (ClinVar) sections are included for
  **context and interest only**. Any medication decision must involve a qualified clinician
  and an accredited clinical-grade test — research-pipeline calls are not validated for care.
- Polygenic scores are reported as a **relative ranking within the requested score set**,
  not as an absolute probability of any outcome. PGS are poorly transferable across
  ancestries and must not be over-interpreted.

## Consent and privacy

- Only analyse genomes you own or for which you have **explicit, informed consent**.
- A genome is identifying and reveals information about biological relatives. The forensic
  STR (CODIS) and kinship modules are genuinely identity-grade — handle outputs accordingly.
- Outputs contain sensitive personal data. Store `--outdir` and `--reference_base` securely
  and do not upload reports to third-party services without consent.

## "Trait" choices

The bundled trait panel (`assets/trait_snps.tsv`) is deliberately limited to non-clinical,
non-stigmatising characteristics (appearance, taste, simple physiology). Behavioural/
psychological SNP associations included are small-effect, much-debated, and provided purely
as conversation pieces — they are **not** deterministic and should not be treated as such.
Adding health-risk variants to this panel is out of scope for the project.
