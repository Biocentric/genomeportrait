process PGS_SCORE_REFERENCE {
    tag "$panel_key"
    label 'process_high'
    // The reference distribution of every score: sample-independent, so computed once.
    storeDir "${params.reference_base}/pgs/${panel_key}/reference_scores"

    conda "bioconda::plink2=2.00a5.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0' :
        'biocontainers/plink2:2.00a5.10--h4ac6f70_0' }"

    input:
    val   panel_key
    tuple path(ref_pgen), path(ref_pvar), path(ref_psam)
    path  scores

    output:
    path "ref_scores",   emit: scores
    path "versions.yml", emit: versions

    script:
    def jobs = Math.max(1, (task.cpus / 2) as int)
    """
    mkdir -p ref_scores logs
    n=0
    for f in ${scores}/*.score; do
        id=\$(basename \$f .score)
        ( plink2 --threads 2 --memory 8000 --pfile ref \\
              --score \$f 1 2 3 header cols=+scoresums --out ref_scores/\$id > logs/\$id.out 2>&1 \\
          || echo "reference scoring failed for \$id" >&2 ) &
        n=\$((n + 1))
        if [ \$((n % ${jobs})) -eq 0 ]; then wait; fi
    done
    wait
    rm -f ref_scores/*.log
    echo "scored: \$(ls ref_scores/*.sscore | wc -l) of \$(ls ${scores}/*.score | wc -l)" >&2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ .*//')
    END_VERSIONS
    """
}
