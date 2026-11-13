process PHARMCAT_REPORT {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::pandas=2.2.2"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? (params.report_sif ?: "${params.reference_base}/containers/genomeportrait-report-${params.report_image_version}.sif") : params.report_docker }"

    input:
    tuple val(meta), path(pharmcat_json)
    path  script   // staged, not on PATH: runs without an exec bit or bin/

    output:
    tuple val(meta), path("*.pgx.tsv"), emit: tsv
    path "versions.yml",                emit: versions

    script:
    """
    python3 ${script} --sample ${meta.id} --json ${pharmcat_json} --out ${meta.id}.pgx.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
