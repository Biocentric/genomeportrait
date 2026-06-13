/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Subworkflow with functionality specific to the your-org/genomeportrait pipeline
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NEXTFLOW_PIPELINE  } from '../../nf-core/utils_nextflow_pipeline'
include { UTILS_NFSCHEMA_PLUGIN    } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap         } from 'plugin/nf-schema'
include { samplesheetToList        } from 'plugin/nf-schema'
include { completionEmail          } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary        } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification           } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE    } from '../../nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:
    ch_versions = Channel.empty()

    UTILS_NEXTFLOW_PIPELINE ( version, true, outdir, workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1 )
    UTILS_NFSCHEMA_PLUGIN   ( workflow, validate_params, null )
    UTILS_NFCORE_PIPELINE   ( nextflow_cli_args )

    //
    // Create channel from input file provided through params.input
    //   Each row -> [ meta(id,sex,platform), fastq_1, fastq_2, bam ]
    //
    Channel
        .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
        .map { meta, fastq_1, fastq_2, bam ->
            // sample/sex/platform are meta columns (see assets/schema_input.json); the
            // positional fields are fastq_1, fastq_2, bam. Backfill meta defaults.
            def m = meta + [ sex: meta.sex ?: 'unknown', platform: meta.platform ?: 'illumina' ]
            return [ m, fastq_1 ?: null, fastq_2 ?: null, bam ?: null ]
        }
        .set { ch_samplesheet }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs, multiqc_reports.getVal())
        }
        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def paramsSummaryMultiqc(summary_params) {
    def summary_section = ''
    summary_params.keySet().each { group ->
        def group_params = summary_params.get(group)
        if (group_params) {
            summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
            summary_section += "    <dl class=\"dl-horizontal\">\n"
            group_params.keySet().sort().each { param ->
                summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999;\">N/A</span>'}</samp></dd>\n"
            }
            summary_section += "    </dl>\n"
        }
    }
    def yaml_file_text = "id: 'genomeportrait-summary'\n"
    yaml_file_text    += "description: ' - this information is collected when the pipeline is started.'\n"
    yaml_file_text    += "section_name: 'genomeportrait Workflow Summary'\n"
    yaml_file_text    += "section_href: 'https://github.com/your-org/genomeportrait'\n"
    yaml_file_text    += "plot_type: 'html'\n"
    yaml_file_text    += "data: |\n"
    yaml_file_text    += "${summary_section}"
    return yaml_file_text
}

def softwareVersionsToYAML(ch_versions) {
    return ch_versions
        .unique()
        .map { version -> version instanceof Path ? version.text : version }
        .unique()
}
