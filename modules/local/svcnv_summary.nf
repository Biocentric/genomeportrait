process SVCNV_SUMMARY {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::pandas=2.2.2"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? params.report_sif : params.report_docker }"

    input:
    tuple val(meta), path(sv_vcf), path(cnr)

    output:
    tuple val(meta), path("*.svcnv_summary.tsv"), emit: tsv
    path "versions.yml",                          emit: versions

    script:
    """
    svcnv_summarise.py \\
        --sample ${meta.id} \\
        ${sv_vcf ? "--sv ${sv_vcf}" : ''} \\
        ${cnr ? "--cnr ${cnr}" : ''} \\
        --out ${meta.id}.svcnv_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
