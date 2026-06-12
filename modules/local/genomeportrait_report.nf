process GENOMEPORTRAIT_REPORT {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::pandas=2.2.2 conda-forge::matplotlib=3.9.2 conda-forge::seaborn=0.13.2 conda-forge::jinja2=3.1.4"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? params.report_sif : params.report_docker }"

    input:
    tuple val(meta), path(parts, stageAs: "parts/*")
    path  multiqc_report

    output:
    tuple val(meta), path("*.genomeportrait.html"), emit: html
    path "versions.yml",                            emit: versions

    script:
    """
    genomeportrait_report.py \\
        --sample ${meta.id} \\
        --parts-dir parts \\
        ${multiqc_report ? "--multiqc ${multiqc_report}" : ''} \\
        --pipeline-version ${workflow.manifest.version} \\
        --out ${meta.id}.genomeportrait.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
