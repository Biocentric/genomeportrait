process STR_REPORT {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::pandas=2.2.2"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? params.report_sif : params.report_docker }"

    input:
    tuple val(meta), path(hipstr_vcf), path(eh_vcf)

    output:
    tuple val(meta), path("${meta.id}_str_profile"), emit: tsv
    path "versions.yml",                         emit: versions

    script:
    """
    str_report.py \\
        --sample ${meta.id} \\
        --hipstr ${hipstr_vcf} \\
        --expansionhunter ${eh_vcf} \\
        --outdir ${meta.id}_str_profile

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
