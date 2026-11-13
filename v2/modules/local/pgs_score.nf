process PGS_SCORE {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::plink2=2.00a5.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0' :
        'biocontainers/plink2:2.00a5.10--h4ac6f70_0' }"

    input:
    tuple val(meta), path(vcf)
    path  afreq
    tuple path(acount), path(allele_wts), path(pca_cols)
    path  scores

    output:
    tuple val(meta), path("target_scores"), path("${meta.id}.pgs_pcs.sscore"), emit: scores
    path "versions.yml",                                                        emit: versions

    script:
    def jobs = Math.max(1, (task.cpus / 2) as int)
    """
    plink2 --threads ${task.cpus} --vcf ${vcf} --autosome --make-pgen --out target

    # Project onto the reference PCs (same command as the reference projection).
    read A1 P1 P2 < ${pca_cols}
    plink2 --threads ${task.cpus} --pfile target --read-freq ${acount} \\
        --score ${allele_wts} 2 \$A1 header-read no-mean-imputation variance-standardize \\
        --score-col-nums \$P1-\$P2 --out ${meta.id}.pgs_pcs

    # Score. --read-freq makes plink2 fill a missing genotype with the reference panel's
    # mean dosage, so an uncovered site contributes exactly what it contributes on average
    # to the reference distribution the sample is compared against.
    mkdir -p target_scores logs
    n=0
    for f in ${scores}/*.score; do
        id=\$(basename \$f .score)
        ( plink2 --threads 2 --memory 4000 --pfile target --read-freq ${afreq} \\
              --score \$f 1 2 3 header cols=+scoresums --out target_scores/\$id > logs/\$id.out 2>&1 \\
          || echo "scoring failed for \$id" >&2 ) &
        n=\$((n + 1))
        if [ \$((n % ${jobs})) -eq 0 ]; then wait; fi
    done
    wait
    rm -f target_scores/*.log target.pgen

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ .*//')
    END_VERSIONS
    """
}
