# your-org/genomeportrait: Output

The pipeline writes everything under `--outdir`. The headline deliverable is
`report/<sample>.genomeportrait.html` — a single self-contained interactive report.

## Directory structure

```
results/
├── report/
│   └── <sample>.genomeportrait.html   # ⭐ the unified portrait (tables + plots)
├── multiqc/
│   └── multiqc_report.html            # technical sequencing/alignment QC
├── qc/                                 # fastp, fastqc, mosdepth, samtools, qualimap
├── alignment/<sample>/                 # markdup + BQSR-recalibrated CRAM
├── variants/
│   ├── deepvariant/<sample>/           # SNV/indel VCF + gVCF
│   ├── manta/<sample>/                 # structural variants
│   └── cnvkit/<sample>/                # copy-number segments
├── annotation/<sample>/                # VEP + gnomAD/ClinVar annotated VCF + TSV
├── ancestry/<sample>/                  # PCA plot, admixture proportions, ROH, kinship
├── lineage/<sample>/                   # mtDNA + Y haplogroups
├── forensic_str/<sample>/              # CODIS STR profile + repeat expansions
├── polygenic/<sample>/                 # per-trait polygenic scores
├── hla/<sample>/                       # HLA class I/II 4-digit type
├── pharmacogenomics/<sample>/          # PharmCAT star alleles + phenotypes
├── traits/<sample>/                    # curated "popular traits" lookup
├── extras/<sample>/                    # telomere length, mtDNA heteroplasmy
└── pipeline_info/                      # execution reports, software versions, params
```

## The GenomePortrait report

`report/<sample>.genomeportrait.html` stitches together every interpretation table and
plot produced for the sample into one navigable document, with a non-medical disclaimer
on every interpretive section. It embeds plots as base64 so the file is fully portable —
email it, archive it, open it offline.

## MultiQC

`multiqc/multiqc_report.html` aggregates read quality (FastQC/fastp), coverage (mosdepth),
alignment stats (samtools/Qualimap) and variant stats (bcftools) across the run.

## Pipeline information

`pipeline_info/` holds the Nextflow execution timeline/report/trace/DAG, the resolved
parameters, and `genomeportrait_software_versions.yml` listing every tool version used.
