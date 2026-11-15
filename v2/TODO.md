# genomeportrait v2 — todo

Status key: `[x]` done and verified by a run or test · `[~]` written, not exercised by a
test · `[ ]` open. See [`PLAN.md`](PLAN.md) for the reasoning and
[`docs/development.md`](docs/development.md) for the runs behind each `[x]`.

## v2.0

### Structure
- [x] Fork `feat/parabricks-pascal` @ `f4eb0ae` into `v2/`
- [x] Samplesheet `gvcf` and `vcf` columns; portrait mode skips alignment and calling (portrait test)
- [x] Remove the redundant `SAMTOOLS_SORT` (smoke test)
- [x] BQSR opt-in; CRAM deliverable without BQSR (smoke test)
- [x] A samplesheet `bam` is called as it is, without MarkDuplicates (modes test)
- [x] References, containers, trimmed FASTQ, markdup BAM and intermediates stay out of `results/`
- [x] Scripts staged as process inputs; no `bin/` on `PATH`
- [x] De-templated (`Biocentric/genomeportrait`, `2.0.0dev`); `nextflow lint` reports no errors
- [x] bwa-mem2 index and BQSR known-sites only prepared when needed (portrait and modes tests)
- [~] `--skip_bqsr false` (BQSR on): unchanged v1 code path, not run in v2

### Polygenic scores
- [x] `assets/pgs/panel_lifestyle.tsv` (204 scores), `assets/pgs/concepts.tsv` (42 concepts)
- [x] PGS Catalog metadata snapshot in `assets/pgs/catalog/`
- [x] Matching, gVCF conversion, scoring, ancestry adjustment (`scripts/pgs_engine.py`)
- [x] Validation against `pgsc_calc`: r = 0.9996 on the same genotypes, r = 0.9997 over the full chain
- [x] gVCF conversion unit test (`tests/test_pgs_target.py`)
- [x] Interpretation reproduces `prs-selfscore` exactly on the same Z values
- [x] Overview page and report fragment
- [x] Nextflow modules and `VCF_POLYGENIC`; cache reused on a second run
- [x] Zero-depth and `--pgs_min_gq` sensitivity measured (`docs/polygenic.md`)
- [~] `--pgs_ids` (extra scores outside the panel): code path written, never run
- [~] `PGS_PANEL_DOWNLOAD`: written, never run — every test used a local copy of the panel

### New analyses
- [x] `GENOME_FACTS`: sex chromosomes, mtDNA copies, depth and breadth (two real genomes, smoke test)
- [~] Resolved sex drives ExpansionHunter and Yleaf: wired and linted; STR and haplogroup
      modules were switched off in both tests

### Report
- [x] Unified report in the overview design; variant overview section now renders
- [x] Standalone `<sample>.pgs_overview.html` and `<sample>.pgs_interpreted.tsv`
- [x] Checked in a browser: light and dark, desktop and 375 px wide, no horizontal overflow

### Tests and docs
- [x] Portrait-mode Nextflow run on a real 21x genome
- [x] `-profile test` smoke test from FASTQ
- [x] All three samplesheet modes in one run on tiny data (`tests/run_test.sh modes`)
- [x] `README.md`, `docs/usage.md`, `docs/output.md`, `docs/polygenic.md`, `docs/ethics.md`, `docs/development.md`, `CHANGELOG.md`
- [ ] Portrait run with **every** interpretation module on (ancestry, lineage, STR, HLA,
      PharmCAT, traits, SV/CNV, VEP). These are v1 modules; v2 changed how nine of them
      receive their script and how two receive the sample's sex. Needs the reference set
      on the RAID, which was not mounted on 2026-11-12.
- [ ] Whole-genome run from FASTQ
- [ ] Commit to a branch and push — nothing in `v2/` is under version control yet

## v2.1

- [x] **ALT-aware alignment** (`--alt_aware`): NA12878 chr22 CYP2D6 locus 2,696 -> 8,935 reads; Cyrius recovers *3/*68+*4
- [ ] Re-align the existing 21x and 42x genomes from FASTQ so CYP2D6 (and HLA/KIR) use ALT-aware alignments
- [ ] Inbreeding F and ROH come from the variant-only VCF (F = -1.0 on the test genome): compute
      them from genotypes that include hom-ref sites (the panel-space target or the ancestry SNPs)
- [ ] KING kinship against the panel sits near -0.05 for an unrelated European; check whether
      mpileup genotyping errors depress it

- [x] **Score the sample in one Python pass** (4 min 50 s → 40 s)
- [x] **Fill missing genotypes from the most similar population's frequencies** (bias up to 0.35 SD removed)
- [ ] Validate the PGS engine on a non-European or admixed genome
- [x] gVCF conversion per chromosome in parallel (5 min 49 s → 1 min 18 s)
- [x] `report_sif`, `hipstr_sif`, `yleaf_sif` follow `reference_base` wherever it is set; cache honours `NXF_SINGULARITY_CACHEDIR`
- [x] Strict opt-in publishing with one `results/<sample>/` tree
- [x] CYP2D6 via Cyrius, passed to PharmCAT as an outside call (6/6 1000 Genomes truth structures)
- [x] Blood groups (ABO, RhD, Kell, Duffy, Kidd), validated on reference-panel frequencies
- [x] One reference panel for ancestry and PGS (hold-out check: unadmixed populations recovered at 0.92–1.00)
- [x] Streamed align → sort → markdup → CRAM (`--markdup samtools`, default)
- [x] VerifyBamID2 contamination estimate (2/5/10% mixtures recovered as 2.2/5.9/10.2%)
- [ ] Stub blocks on every process; nf-test

## v2.2
- [x] 23andMe-format export (99.99% agreement with the VCF; plink --23file reads it)
- [x] Archaic introgression share (hmmix; 8 1000 Genomes genomes vs the published calls: 87-92% recovered, 91-92% precise, labels 99-100% agree)
- [x] Chromosome painting (Beagle + RFMix; 20 hold-out 1000 Genomes genomes: unadmixed 97-100% own group, ASW 77/21% and ACB 88/11% African/European as published)
- [x] KIR typing (T1K; 59 1000 Genomes genomes: framework genes in all, B-haplotype gene frequencies as published)
