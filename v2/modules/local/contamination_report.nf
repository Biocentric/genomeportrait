process CONTAMINATION_REPORT {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? (params.report_sif ?: "${params.reference_base}/containers/genomeportrait-report-${params.report_image_version}.sif") : params.report_docker }"

    input:
    tuple val(meta), path(selfsm)
    path  script

    output:
    tuple val(meta), path("${meta.id}.contamination_facts.tsv"), emit: tsv
    path "versions.yml",                                         emit: versions

    script:
    """
    python3 ${script} --sample ${meta.id} --selfsm ${selfsm} --out ${meta.id}.contamination_facts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
