process ANCESTRY_REPORT {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::pandas=2.2.2 conda-forge::matplotlib=3.9.2"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? params.report_sif : params.report_docker }"

    input:
    tuple val(meta), path(eigenvec), path(admixture_q), path(roh), path(kin)
    tuple path(ref_pgen), path(ref_pvar), path(ref_psam)

    output:
    tuple val(meta), path("${meta.id}_ancestry"), emit: bundle
    path "versions.yml",                          emit: versions

    script:
    """
    ancestry_report.py \\
        --sample ${meta.id} \\
        --eigenvec ${eigenvec} \\
        --admixture ${admixture_q} \\
        --roh ${roh} \\
        --kin ${kin} \\
        --ref-psam ${ref_psam} \\
        --outdir ${meta.id}_ancestry

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
