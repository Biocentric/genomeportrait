# genomeportrait: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - [unreleased]

Initial release of genomeportrait, created following the nf-core template.

### `Added`

- Reference auto-download & caching subworkflow (GRCh38, VEP cache, gnomAD v4, ClinVar, 1000G+HGDP panel, PGS Catalog, STR catalogs, PharmCAT, HLA resources).
- Read QC/trimming (fastp, FastQC) and alignment (bwa-mem2 / minimap2).
- Duplicate marking + BQSR (GATK4).
- Alignment QC (mosdepth, samtools stats, Qualimap).
- Small-variant calling with DeepVariant + GLnexus.
- SV/CNV calling with Manta + CNVkit.
- Annotation with Ensembl VEP + SnpSift (gnomAD/ClinVar/dbNSFP).
- Ancestry/admixture (PLINK2 + ADMIXTURE + PCA), ROH/inbreeding, KING kinship.
- mtDNA (HaploGrep2) and Y (Yleaf) haplogroups.
- Forensic STR profile (HipSTR CODIS + ExpansionHunter).
- Polygenic scoring (pgsc_calc-style PLINK2 --score, PGS Catalog).
- HLA typing (T1K).
- Pharmacogenomics (PharmCAT).
- Curated non-clinical "popular traits" SNP lookup.
- Telomere length (Telseq) and mtDNA heteroplasmy (mutserve).
- Unified GenomePortrait HTML report + MultiQC.

### `Fixed`

### `Dependencies`

### `Deprecated`
