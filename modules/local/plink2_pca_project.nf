process PLINK2_PCA_PROJECT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::plink2=2.00a5.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0' :
        'biocontainers/plink2:2.00a5.10--h4ac6f70_0' }"

    input:
    tuple val(meta), path(pgen), path(pvar), path(psam)
    tuple path(ref_pgen), path(ref_pvar), path(ref_psam)

    output:
    tuple val(meta), path("*.eigenvec"),      emit: eigenvec
    tuple val(meta), path("*.eigenvec.allele"), emit: loadings, optional: true
    path "versions.yml",                      emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    # Learn PCA on the reference panel, then project the merged sample onto it
    plink2 --pfile ${ref_pgen.baseName} $args \\
        --pca 10 allele-wts --out ref_pca
    plink2 --pfile ${pgen.baseName} \\
        --read-freq ref_pca.afreq \\
        --score ref_pca.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \\
        --score-col-nums 6-15 --out ${meta.id}.proj
    # Reformat scores as an eigenvec table (PC1..PC10 per individual)
    cp ${meta.id}.proj.sscore ${meta.id}.eigenvec

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ .*//')
    END_VERSIONS
    """
}
