process EXTRAS_REPORT {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::pandas=2.2.2"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? (params.report_sif ?: "${params.reference_base}/containers/genomeportrait-report-${params.report_image_version}.sif") : params.report_docker }"

    input:
    tuple val(meta), path(telseq), path(mito)
    path  script   // staged, not on PATH: runs without an exec bit or bin/

    output:
    tuple val(meta), path("*.extras.tsv"), emit: tsv
    path "versions.yml",                   emit: versions

    script:
    """
    python3 ${script} \\
        --sample ${meta.id} \\
        ${telseq ? "--telseq ${telseq}" : ''} \\
        ${mito ? "--mito ${mito}" : ''} \\
        --out ${meta.id}.extras.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.id}.extras.tsv"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: true
    END_VERSIONS
    """
}
