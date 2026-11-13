process PGS_SITES {
    tag "$panel_key"
    label 'process_medium'
    storeDir "${params.reference_base}/pgs/${panel_key}/sites"

    conda "conda-forge::python=3.12 conda-forge::numpy=2.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? params.report_sif : params.report_docker }"

    input:
    val   panel_key
    path  scorefiles
    path  engine

    output:
    path "pgs_sites.bed1", emit: sites
    path "versions.yml",   emit: versions

    script:
    """
    python3 ${engine} sites --scorefiles ${scorefiles} --out pgs_sites.bed1 --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
