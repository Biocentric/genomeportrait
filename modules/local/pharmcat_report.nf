process PHARMCAT_REPORT {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::pandas=2.2.2"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? params.report_sif : params.report_docker }"

    input:
    tuple val(meta), path(pharmcat_json)

    output:
    tuple val(meta), path("*.pgx.tsv"), emit: tsv
    path "versions.yml",                emit: versions

    script:
    """
    pharmcat_report.py --sample ${meta.id} --json ${pharmcat_json} --out ${meta.id}.pgx.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
