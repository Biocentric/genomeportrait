process PLINK2_VCF2PGEN {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::plink2=2.00a5.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0' :
        'biocontainers/plink2:2.00a5.10--h4ac6f70_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.id}.pgen"), path("${meta.id}.pvar"), path("${meta.id}.psam"), emit: pgen
    path "versions.yml", emit: versions

    script:
    """
    # ID = chr:pos so PGS scoring files (matched on position + effect allele) can align;
    # drop duplicate positions (split multiallelics) keeping the first.
    plink2 --vcf $vcf --set-all-var-ids '@:#' --rm-dup force-first \\
        --make-pgen --out ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ .*//')
    END_VERSIONS
    """
}
