process TRAIT_LOOKUP {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::pandas=2.2.2"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? params.report_sif : params.report_docker }"

    input:
    tuple val(meta), path(geno)
    path  panel

    output:
    tuple val(meta), path("*.traits.tsv"), emit: tsv
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    trait_lookup.py \\
        --sample ${meta.id} \\
        --geno ${geno} \\
        --panel ${panel} \\
        --out ${meta.id}.traits.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
