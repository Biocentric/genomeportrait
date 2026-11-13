process PGS_MATCH {
    tag "$panel_key"
    label 'process_medium'
    storeDir "${params.reference_base}/pgs/${panel_key}/match"

    conda "conda-forge::python=3.12 conda-forge::numpy=2.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? (params.report_sif ?: "${params.reference_base}/containers/genomeportrait-report-${params.report_image_version}.sif") : params.report_docker }"

    input:
    val   panel_key
    path  scorefiles
    tuple path(ref_pgen), path(ref_pvar), path(ref_psam)
    path  afreq
    path  engine

    output:
    path "scores",              emit: scores
    path "pgs_match_stats.tsv", emit: stats
    path "versions.yml",        emit: versions

    script:
    // Matching is by chr:pos AND alleles. A position-only check once reported 99.95%
    // overlap for a target in which not a single allele could match.
    def amb = params.pgs_keep_ambiguous ? '--keep-ambiguous' : ''
    """
    python3 ${engine} match --scorefiles ${scorefiles} --pvar ${ref_pvar} --afreq ${afreq} \\
        --outdir scores --stats pgs_match_stats.tsv ${amb} --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
