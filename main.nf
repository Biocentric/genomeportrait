#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    your-org/genomeportrait
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    State-of-the-art personal whole-genome sequencing analysis.
    Github : https://github.com/your-org/genomeportrait
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GENOMEPORTRAIT      } from './workflows/genomeportrait'
include { PREPARE_GENOME      } from './subworkflows/local/prepare_genome'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_genomeportrait_pipeline'
include { PIPELINE_COMPLETION      } from './subworkflows/local/utils_genomeportrait_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow YOURORG_GENOMEPORTRAIT {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // SUBWORKFLOW: Download and index all reference resources (cached) + build report container
    //
    PREPARE_GENOME ( params.genome )

    // Assemble a plain map of value-channels. Property access (ch_reference.fasta)
    // therefore returns a reusable value-channel with proper data dependencies.
    ch_reference = [
        fasta         : PREPARE_GENOME.out.fasta,
        fai           : PREPARE_GENOME.out.fai,
        dict          : PREPARE_GENOME.out.dict,
        bwamem2_index : PREPARE_GENOME.out.bwamem2_index,
        dbsnp         : PREPARE_GENOME.out.dbsnp,
        known_indels  : PREPARE_GENOME.out.known_indels,
        mills         : PREPARE_GENOME.out.mills,
        vep_cache     : PREPARE_GENOME.out.vep_cache,
        gnomad        : PREPARE_GENOME.out.gnomad,
        clinvar       : PREPARE_GENOME.out.clinvar,
        kgp_hgdp_panel: PREPARE_GENOME.out.kgp_hgdp,
        str_catalog   : PREPARE_GENOME.out.str_catalog,
        hipstr_codis  : PREPARE_GENOME.out.hipstr_codis,
        report_sif    : PREPARE_GENOME.out.report_sif,
        hipstr_ready  : PREPARE_GENOME.out.hipstr_ready,
        hla_index     : PREPARE_GENOME.out.hla_index,
        yleaf_data    : PREPARE_GENOME.out.yleaf_data,
        yleaf_ready   : PREPARE_GENOME.out.yleaf_ready
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
    YOURORG_GENOMEPORTRAIT ( PIPELINE_INITIALISATION.out.samplesheet )

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
        YOURORG_GENOMEPORTRAIT.out.multiqc_report
    )
}
