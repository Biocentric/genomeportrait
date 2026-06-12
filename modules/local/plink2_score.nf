process PLINK2_SCORE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::plink2=2.00a5.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0' :
        'biocontainers/plink2:2.00a5.10--h4ac6f70_0' }"

    input:
    tuple val(meta), path(pgen), path(pvar), path(psam)
    path  scorefile

    output:
    tuple val(meta), path("*.sscore"), emit: sscore
    path "versions.yml",               emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    # scorefile columns: 1=variant_id (chr:pos:ref:alt), 2=effect_allele, 3..=per-PGS weights
    zcat -f $scorefile > scorefile.txt
    NCOL=\$(head -n1 scorefile.txt | awk '{print NF}')
    plink2 --pfile ${pgen.baseName} \\
        --score scorefile.txt 1 2 header-read cols=+scoresums,+denom ignore-dup-ids \\
        --score-col-nums 3-\${NCOL} \\
        --out ${meta.id}.prs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ .*//')
    END_VERSIONS
    """
}
