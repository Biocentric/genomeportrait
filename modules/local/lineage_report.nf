process LINEAGE_REPORT {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::pandas=2.2.2"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? params.report_sif : params.report_docker }"

    input:
    tuple val(meta), path(mt_result), path(y_result)

    output:
    tuple val(meta), path("*.lineage.tsv"), emit: tsv
    path "versions.yml",                    emit: versions

    script:
    """
    lineage_report.py \\
        --sample ${meta.id} \\
        --mtdna ${mt_result} \\
        ${y_result ? "--ydna ${y_result}" : ''} \\
        --out ${meta.id}.lineage.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
