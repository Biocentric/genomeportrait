#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Biocentric/genomeportrait
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Personal whole-genome portrait, version 2.
    Github : https://github.com/Biocentric/genomeportrait
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GENOMEPORTRAIT      } from './workflows/genomeportrait'
include { samplesheetToList   } from 'plugin/nf-schema'
include { PREPARE_GENOME      } from './subworkflows/local/prepare_genome'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_genomeportrait_pipeline'
include { PIPELINE_COMPLETION      } from './subworkflows/local/utils_genomeportrait_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BIOCENTRIC_GENOMEPORTRAIT {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // SUBWORKFLOW: Download and index all reference resources (cached) + build report container
    //
    // Only build what the samplesheet needs: a run that starts from alignments or from
    // finished call sets has no use for the bwa-mem2 index (~17 GB, ~90 min to build).
    def rows = samplesheetToList(params.input, "${projectDir}/assets/schema_input.json")
    def needs_alignment = rows.any { r -> r[3] == null || r[3] == [] }
    PREPARE_GENOME ( params.genome, needs_alignment )

    // Assemble a plain map of value-channels. Property access (ch_reference.fasta)
    // therefore returns a reusable value-channel with proper data dependencies.
    ch_reference = [
        fasta         : PREPARE_GENOME.out.fasta,
        fai           : PREPARE_GENOME.out.fai,
        dict          : PREPARE_GENOME.out.dict,
        bwamem2_index : PREPARE_GENOME.out.bwamem2_index,
        bwa_alt       : PREPARE_GENOME.out.bwa_alt,
        dbsnp         : PREPARE_GENOME.out.dbsnp,
        known_indels  : PREPARE_GENOME.out.known_indels,
        mills         : PREPARE_GENOME.out.mills,
        vep_cache     : PREPARE_GENOME.out.vep_cache,
        gnomad        : PREPARE_GENOME.out.gnomad,
        clinvar       : PREPARE_GENOME.out.clinvar,
        panel         : PREPARE_GENOME.out.panel,         // HGDP+1kGP: ancestry and PGS
        str_catalog   : PREPARE_GENOME.out.str_catalog,
        hipstr_codis  : PREPARE_GENOME.out.hipstr_codis,
        report_sif    : PREPARE_GENOME.out.report_sif,
        hipstr_ready  : PREPARE_GENOME.out.hipstr_ready,
        hla_index     : PREPARE_GENOME.out.hla_index,
        kir_index     : PREPARE_GENOME.out.kir_index,
        yleaf_data    : PREPARE_GENOME.out.yleaf_data,
        yleaf_ready   : PREPARE_GENOME.out.yleaf_ready,
        cyrius_ready  : PREPARE_GENOME.out.cyrius_ready
    ]

    //
    // WORKFLOW: Run the main analysis pipeline
    //
    GENOMEPORTRAIT (
        samplesheet,
        ch_reference
    )

    emit:
    multiqc_report = GENOMEPORTRAIT.out.multiqc_report
    portrait       = GENOMEPORTRAIT.out.portrait
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    BIOCENTRIC_GENOMEPORTRAIT ( PIPELINE_INITIALISATION.out.samplesheet )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        BIOCENTRIC_GENOMEPORTRAIT.out.multiqc_report
    )
}
