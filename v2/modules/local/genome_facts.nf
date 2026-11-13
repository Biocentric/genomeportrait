process GENOME_FACTS {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? (params.report_sif ?: "${params.reference_base}/containers/genomeportrait-report-${params.report_image_version}.sif") : params.report_docker }"

    input:
    tuple val(meta), path(summary), path(dist)
    path  script

    output:
    tuple val(meta), path("${meta.id}.genome_facts.tsv"), emit: tsv
    tuple val(meta), path("${meta.id}.sex.txt"),          emit: sex
    path "versions.yml",                                  emit: versions

    script:
    """
    python3 ${script} --sample ${meta.id} --summary ${summary} --dist ${dist} \\
        --declared ${meta.sex ?: 'unknown'} \\
        --out ${meta.id}.genome_facts.tsv --sex-out ${meta.id}.sex.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
