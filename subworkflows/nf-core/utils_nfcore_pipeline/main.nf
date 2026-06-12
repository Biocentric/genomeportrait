//
// Subworkflow with utility functions specific to the nf-core pipeline template
//   (compact local copy of the nf-core/utils_nfcore_pipeline subworkflow)
//

workflow UTILS_NFCORE_PIPELINE {
    take:
    nextflow_cli_args

    main:
    valid_config = checkConfigProvided()

    emit:
    valid_config
}

def checkConfigProvided() {
    def valid_config = true as Boolean
    if (workflow.profile == 'standard' && workflow.configFiles.size() <= 1) {
        log.warn "[${workflow.manifest.name}] You are attempting to run the pipeline without any custom configuration!\n" +
            "This will be dependent on your local compute environment but can be achieved via one or more of the following:\n" +
            "   (1) Using an existing pipeline shared profile (e.g. -profile docker or -profile singularity)\n" +
            "   (2) Using an existing nf-core/configs institutional profile\n" +
            "   (3) Using your own local custom config (e.g. -c run.config)"
        valid_config = false
    }
    return valid_config
}

//
// Generate version string
//
def getWorkflowVersion() {
    def version_string = "" as String
    if (workflow.manifest.version) {
        def prefix_v = workflow.manifest.version[0] != 'v' ? 'v' : ''
        version_string += "${prefix_v}${workflow.manifest.version}"
    }
    if (workflow.commitId) {
        version_string += "-g${workflow.commitId.substring(0, 7)}"
    }
    return version_string
}

//
// Get software versions for the pipeline
//
def processVersionsFromYAML(yaml_file) {
    def yaml = new org.yaml.snakeyaml.Yaml()
    def versions = yaml.load(yaml_file).collectEntries { k, v -> [k.tokenize(':')[-1], v] }
    return yaml.dumpAsMap(versions).trim()
}

def workflowVersionToYAML() {
    return """
    Workflow:
        ${workflow.manifest.name}: ${getWorkflowVersion()}
        Nextflow: ${workflow.nextflow.version}
    """.stripIndent().trim()
}

def softwareVersionsToYAML(ch_versions) {
    return ch_versions
        .unique()
        .map { version -> processVersionsFromYAML(version) }
        .unique()
        .mix(Channel.of(workflowVersionToYAML()))
}

//
// Completion email + summary (compact)
//
def completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs, multiqc_report) {
    def subject = "[${workflow.manifest.name}] ${workflow.success ? 'Successful' : 'FAILED'}: ${workflow.runName}"
    def email_address = workflow.success ? email : (email_on_fail ?: email)
    if (!email_address) return
    log.info "[${workflow.manifest.name}] Sending completion email to ${email_address}"
    def msg = "Pipeline ${workflow.manifest.name} run ${workflow.runName} ${workflow.success ? 'completed successfully' : 'FAILED'}.\nDuration: ${workflow.duration}\nResults: ${outdir}"
    try {
        ['mail', '-s', subject, email_address].execute() << msg
    } catch (all) {
        log.warn "[${workflow.manifest.name}] Could not send mail (is `mail` installed?). Summary logged instead."
    }
}

def completionSummary(monochrome_logs) {
    if (workflow.success) {
        log.info "-\033[0;32m[${workflow.manifest.name}] Pipeline completed successfully\033[0m-"
    } else {
        log.info "-\033[0;31m[${workflow.manifest.name}] Pipeline completed with errors\033[0m-"
    }
}

def imNotification(summary_params, hook_url) {
    def msg = [text: "[${workflow.manifest.name}] ${workflow.runName} ${workflow.success ? 'succeeded' : 'failed'}"]
    try {
        def http = new URL(hook_url).openConnection() as java.net.HttpURLConnection
        http.setRequestMethod('POST')
        http.setRequestProperty('Content-Type', 'application/json')
        http.setDoOutput(true)
        http.outputStream.withWriter { it << groovy.json.JsonOutput.toJson(msg) }
        http.responseCode
    } catch (all) {
        log.warn "[${workflow.manifest.name}] Could not POST notification to hook_url."
    }
}
