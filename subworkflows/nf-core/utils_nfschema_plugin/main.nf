//
// Subworkflow that uses the nf-schema plugin to validate parameters and render help
//   (compact local copy of the nf-core/utils_nfschema_plugin subworkflow)
//

include { paramsSummaryLog   } from 'plugin/nf-schema'
include { validateParameters } from 'plugin/nf-schema'

workflow UTILS_NFSCHEMA_PLUGIN {
    take:
    input_workflow      // the workflow object
    validate_params     // boolean: validate the parameters
    parameters_schema   // string : path to the parameters JSON schema

    main:
    def schema_filename = parameters_schema ?: "nextflow_schema.json"
    log.info paramsSummaryLog(input_workflow, parameters_schema: schema_filename)
    if (validate_params) {
        validateParameters(parameters_schema: schema_filename)
    }

    emit:
    dummy_emit = true
}
