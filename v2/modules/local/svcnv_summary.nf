process SVCNV_SUMMARY {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::pandas=2.2.2"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? (params.report_sif ?: "${params.reference_base}/containers/genomeportrait-report-${params.report_image_version}.sif") : params.report_docker }"

    input:
    tuple val(meta), path(sv_vcf), path(cnr)
    path  script   // staged, not on PATH: runs without an exec bit or bin/

    output:
    tuple val(meta), path("${meta.id}_svcnv_summary"), emit: tsv
    path "versions.yml",                          emit: versions

    script:
    """
    python3 ${script} \\
        --sample ${meta.id} \\
        ${sv_vcf ? "--sv ${sv_vcf}" : ''} \\
        ${cnr ? "--cnr ${cnr}" : ''} \\
        --outdir ${meta.id}_svcnv_summary

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p "${meta.id}_svcnv_summary"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: true
    END_VERSIONS
    """
}
