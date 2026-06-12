process EXTRAS_REPORT {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::pandas=2.2.2"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? params.report_sif : params.report_docker }"

    input:
    tuple val(meta), path(telseq), path(mito)

    output:
    tuple val(meta), path("*.extras.tsv"), emit: tsv
    path "versions.yml",                   emit: versions

    script:
    """
    extras_report.py \\
        --sample ${meta.id} \\
        ${telseq ? "--telseq ${telseq}" : ''} \\
        ${mito ? "--mito ${mito}" : ''} \\
        --out ${meta.id}.extras.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
