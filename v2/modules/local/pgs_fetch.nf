process PGS_FETCH {
    tag "$panel_key"
    label 'process_low'
    // Scoring files and catalogue metadata depend only on the panel, never on the sample.
    storeDir "${params.reference_base}/pgs/${panel_key}/fetch"

    conda "conda-forge::python=3.12 conda-forge::numpy=2.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? params.report_sif : params.report_docker }"

    input:
    val   panel_key
    path  panel
    val   extra_ids
    path  catalog_seed, stageAs: 'catalog_seed'   // never the same name as the 'catalog' output
    path  fetch_script
    path  catalog_script

    output:
    path "scorefiles",   emit: scorefiles
    path "catalog",      emit: catalog
    path "versions.yml", emit: versions

    script:
    def extra = extra_ids ? "--extra ${extra_ids}" : ''
    def ids   = extra_ids ? "--ids ${extra_ids}"   : ''
    """
    python3 ${fetch_script} --panel ${panel} ${extra} --outdir scorefiles --jobs 4
    # Metadata for panel scores ships with the pipeline; only scores added with --pgs_ids
    # are looked up, so a default run needs no access to the PGS Catalog API.
    python3 ${catalog_script} --panel ${panel} ${ids} --seed catalog_seed --outdir catalog

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
