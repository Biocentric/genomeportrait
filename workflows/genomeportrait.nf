/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// QC / prep
include { FASTQ_QC_TRIM            } from '../subworkflows/local/fastq_qc_trim'
include { CHECK_GPU                } from '../modules/local/check_gpu'
include { ALIGN_AND_CALL           } from '../subworkflows/local/align_and_call'
include { BAM_QC                   } from '../subworkflows/local/bam_qc'

// Variants
include { BAM_CALL_SV_CNV          } from '../subworkflows/local/bam_call_sv_cnv'
include { VCF_ANNOTATE             } from '../subworkflows/local/vcf_annotate'

// Interpretation layers
include { BAM_ANCESTRY             } from '../subworkflows/local/vcf_ancestry'
include { BAM_HAPLOGROUPS          } from '../subworkflows/local/bam_haplogroups'
include { BAM_FORENSIC_STR         } from '../subworkflows/local/bam_forensic_str'
include { VCF_POLYGENIC            } from '../subworkflows/local/vcf_polygenic'
include { BAM_HLA                  } from '../subworkflows/local/bam_hla'
include { VCF_PHARMACOGENOMICS     } from '../subworkflows/local/vcf_pharmacogenomics'
include { VCF_TRAITS               } from '../subworkflows/local/vcf_traits'
include { BAM_EXTRAS               } from '../subworkflows/local/bam_extras'

// Reporting
include { GENOMEPORTRAIT_REPORT    } from '../modules/local/genomeportrait_report'
include { MULTIQC                  } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap         } from 'plugin/nf-schema'
include { paramsSummaryMultiqc     } from '../subworkflows/local/utils_genomeportrait_pipeline'
include { softwareVersionsToYAML   } from '../subworkflows/local/utils_genomeportrait_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENOMEPORTRAIT {

    take:
    ch_samplesheet // channel: [ meta, fastq_1, fastq_2, bam ]
    ch_reference   // map    : prepared reference resources from PREPARE_GENOME

    main:

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()
    // Per-sample result collectors feeding the final portrait report
    ch_report_parts  = Channel.empty()

    //
    // Branch inputs: start from FASTQ or from existing alignment
    //
    ch_samplesheet
        .branch { meta, fq1, fq2, bam ->
            fastq: bam == null || bam == []
                return [ meta, [ fq1, fq2 ].findAll { it } ]
            aligned: true
                return [ meta, bam ]
        }
        .set { ch_inputs }

    //
    // STAGE 1: Read QC + trimming
    //
    FASTQ_QC_TRIM ( ch_inputs.fastq, params.skip_trimming )
    ch_versions      = ch_versions.mix(FASTQ_QC_TRIM.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_QC_TRIM.out.multiqc.collect{it[1]})

    //
    // STAGE 1b: Decide GPU vs CPU toolchain (Parabricks when a capable GPU is present)
    //
    CHECK_GPU ()
    ch_gpu_ok   = CHECK_GPU.out.decision.map { it.text.trim() == 'true' }.first()
    ch_versions = ch_versions.mix(CHECK_GPU.out.versions)
    CHECK_GPU.out.report.subscribe { log.info "[genomeportrait] GPU check:\n${it.text.trim()}" }

    //
    // STAGE 2: Align + mark-dup/BQSR + small-variant calling (GPU Parabricks or CPU path)
    //
    ch_prealigned = ch_inputs.aligned.map { meta, bam -> [ meta, bam, "${bam}.bai" ] }
    ALIGN_AND_CALL ( FASTQ_QC_TRIM.out.reads, ch_prealigned, ch_reference, ch_gpu_ok )

    // Gate all downstream stages on the report container being built (avoids resume race)
    ch_build_ready   = ch_reference.report_sif.ifEmpty('ready').first()
    ch_analysis_bam  = ALIGN_AND_CALL.out.bam
        .combine( ch_build_ready )
        .map { meta, bam, bai, ready -> [ meta, bam, bai ] }
    ch_vcf           = ALIGN_AND_CALL.out.vcf      // [ meta, vcf.gz, tbi ]
    ch_versions      = ch_versions.mix(ALIGN_AND_CALL.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(ALIGN_AND_CALL.out.metrics.collect{it[1]})
    ch_report_parts  = ch_report_parts.mix(ALIGN_AND_CALL.out.stats.map{ m,f -> [m,'variant_stats',f] })

    //
    // STAGE 3: Alignment QC
    //
    if (!params.skip_qc) {
        BAM_QC ( ch_analysis_bam, BAM_MARKDUP_BQSR.out.md_bam, ch_reference )
        ch_versions      = ch_versions.mix(BAM_QC.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_QC.out.multiqc.collect{it[1]})
    }

    //
    // STAGE 5: SV / CNV
    //
    if (params.call_sv || params.call_cnv) {
        BAM_CALL_SV_CNV ( ch_analysis_bam, ch_reference )
        ch_versions     = ch_versions.mix(BAM_CALL_SV_CNV.out.versions)
        ch_report_parts = ch_report_parts.mix(BAM_CALL_SV_CNV.out.summary.map{ m,f -> [m,'sv_cnv',f] })
    }

    //
    // STAGE 6: Annotation (VEP + gnomAD/ClinVar/dbNSFP)
    //
    if (!params.skip_annotation) {
        VCF_ANNOTATE ( ch_vcf, ch_reference )
        ch_vcf_annotated = VCF_ANNOTATE.out.vcf
        ch_versions      = ch_versions.mix(VCF_ANNOTATE.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(VCF_ANNOTATE.out.reports.collect{it[1]})
        ch_report_parts  = ch_report_parts.mix(VCF_ANNOTATE.out.tsv.map{ m,f -> [m,'annotation',f] })
    } else {
        ch_vcf_annotated = ch_vcf
    }

    //
    // STAGE 7: Ancestry & admixture (opt-in; needs 1000G+HGDP panel)
    //
    if (!params.skip_ancestry) {
        BAM_ANCESTRY ( ch_analysis_bam, ch_vcf, ch_reference )
        ch_versions     = ch_versions.mix(BAM_ANCESTRY.out.versions)
        ch_report_parts = ch_report_parts.mix(BAM_ANCESTRY.out.results.map{ m,f -> [m,'ancestry',f] })
    }

    //
    // STAGE 8: Maternal / paternal lineage haplogroups
    //
    if (!params.skip_haplogroups) {
        BAM_HAPLOGROUPS ( ch_analysis_bam, ch_vcf, ch_reference )
        ch_versions     = ch_versions.mix(BAM_HAPLOGROUPS.out.versions)
        ch_report_parts = ch_report_parts.mix(BAM_HAPLOGROUPS.out.results.map{ m,f -> [m,'haplogroups',f] })
    }

    //
    // STAGE 9: Forensic STR profile (CODIS) + repeat expansions
    //
    if (!params.skip_str) {
        BAM_FORENSIC_STR ( ch_analysis_bam, ch_reference )
        ch_versions     = ch_versions.mix(BAM_FORENSIC_STR.out.versions)
        ch_report_parts = ch_report_parts.mix(BAM_FORENSIC_STR.out.results.map{ m,f -> [m,'forensic_str',f] })
    }

    //
    // STAGE 10: Polygenic scores (PGS Catalog)
    //
    if (!params.skip_prs) {
        VCF_POLYGENIC ( ch_vcf, ch_reference )
        ch_versions     = ch_versions.mix(VCF_POLYGENIC.out.versions)
        ch_report_parts = ch_report_parts.mix(VCF_POLYGENIC.out.scores.map{ m,f -> [m,'polygenic',f] })
    }

    //
    // STAGE 11: HLA typing (T1K)
    //
    if (!params.skip_hla) {
        BAM_HLA ( ch_analysis_bam, ch_reference )
        ch_versions     = ch_versions.mix(BAM_HLA.out.versions)
        ch_report_parts = ch_report_parts.mix(BAM_HLA.out.results.map{ m,f -> [m,'hla',f] })
    }

    //
    // STAGE 12: Pharmacogenomics (PharmCAT)
    //
    if (!params.skip_pharmcat) {
        VCF_PHARMACOGENOMICS ( ch_vcf, ch_reference )
        ch_versions     = ch_versions.mix(VCF_PHARMACOGENOMICS.out.versions)
        ch_report_parts = ch_report_parts.mix(VCF_PHARMACOGENOMICS.out.results.map{ m,f -> [m,'pharmcat',f] })
    }

    //
    // STAGE 13: "Popular" non-clinical trait lookup
    //
    if (!params.skip_traits) {
        VCF_TRAITS ( ch_analysis_bam, ch_reference, file(params.trait_panel) )
        ch_versions     = ch_versions.mix(VCF_TRAITS.out.versions)
        ch_report_parts = ch_report_parts.mix(VCF_TRAITS.out.results.map{ m,f -> [m,'traits',f] })
    }

    //
    // STAGE 14: Extras (telomere length, mtDNA heteroplasmy)
    //
    BAM_EXTRAS ( ch_analysis_bam, ch_reference )
    ch_versions     = ch_versions.mix(BAM_EXTRAS.out.versions)
    ch_report_parts = ch_report_parts.mix(BAM_EXTRAS.out.results.map{ m,f -> [m,'extras',f] })

    //
    // Collate software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'genomeportrait_software_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // STAGE 15a: MultiQC
    //
    ch_multiqc_report = Channel.empty()
    if (!params.skip_multiqc) {
        ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
        summary_params    = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
        ch_multiqc_files = ch_multiqc_files
            .mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
            .mix(ch_collated_versions)

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            [], [], [], []
        )
        ch_multiqc_report = MULTIQC.out.report
        ch_versions       = ch_versions.mix(MULTIQC.out.versions)
    }

    //
    // STAGE 15b: Unified GenomePortrait HTML report
    //   Group all per-sample report parts into one bundle per sample.
    //
    ch_report_parts
        .map { meta, tag, f -> [ meta, [ tag, f ] ] }
        .groupTuple()
        .map { meta, parts -> [ meta, parts.collect{ it[1] } ] }
        .set { ch_report_bundle }

    GENOMEPORTRAIT_REPORT (
        ch_report_bundle,
        ch_multiqc_report.ifEmpty([])
    )
    ch_versions = ch_versions.mix(GENOMEPORTRAIT_REPORT.out.versions)

    emit:
    multiqc_report = ch_multiqc_report
    portrait       = GENOMEPORTRAIT_REPORT.out.html
    versions       = ch_versions
}
