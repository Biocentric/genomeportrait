process TRAIT_LOOKUP {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::pandas=2.2.2"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? (params.report_sif ?: "${params.reference_base}/containers/genomeportrait-report-${params.report_image_version}.sif") : params.report_docker }"

    input:
    tuple val(meta), path(geno)
    path  panel
    path  script   // staged, not on PATH: runs without an exec bit or bin/

    output:
    tuple val(meta), path("*.traits.tsv"), emit: tsv
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python3 ${script} \\
        --sample ${meta.id} \\
        --geno ${geno} \\
        --panel ${panel} \\
        --out ${meta.id}.traits.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.id}.traits.tsv"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: true
    END_VERSIONS
    """
}
