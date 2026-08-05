process PRS_REPORT {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::pandas=2.2.2 conda-forge::scipy=1.14.1 conda-forge::matplotlib=3.9.2"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? params.report_sif : params.report_docker }"

    input:
    tuple val(meta), path(sscore)
    path  metadata
    path  ref_scores

    output:
    tuple val(meta), path("*.prs.tsv"), emit: tsv
    path "versions.yml",                emit: versions

    script:
    """
    prs_report.py \\
        --sample ${meta.id} \\
        --indir . \\
        --metadata ${metadata} \\
        ${ref_scores ? "--ref-scores ${ref_scores}" : ''} \\
        --min-overlap ${params.pgs_min_overlap} \\
        --out ${meta.id}.prs.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
