process HLA_REPORT {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::pandas=2.2.2"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? params.report_sif : params.report_docker }"

    input:
    tuple val(meta), path(t1k_genotype)

    output:
    tuple val(meta), path("*.hla.tsv"), emit: tsv
    path "versions.yml",                emit: versions

    script:
    """
    hla_report.py --sample ${meta.id} --t1k ${t1k_genotype} --out ${meta.id}.hla.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
