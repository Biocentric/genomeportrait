process PGS_PROJECT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::plink2=2.00a5.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0' :
        'biocontainers/plink2:2.00a5.10--h4ac6f70_0' }"

    input:
    tuple val(meta), path(vcf)
    tuple path(acount), path(allele_wts), path(pca_cols)

    output:
    tuple val(meta), path("${meta.id}.pgs_pcs.sscore"), emit: pcs
    path "versions.yml",                                emit: versions

    script:
    """
    plink2 --threads ${task.cpus} --vcf ${vcf} --autosome --make-pgen --out target

    # Project onto the reference PCs with the same command that projected the reference.
    read A1 P1 P2 < ${pca_cols}
    plink2 --threads ${task.cpus} --pfile target --read-freq ${acount} \\
        --score ${allele_wts} 2 \$A1 header-read no-mean-imputation variance-standardize \\
        --score-col-nums \$P1-\$P2 --out ${meta.id}.pgs_pcs
    rm -f target.pgen

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ .*//')
    END_VERSIONS
    """
}
