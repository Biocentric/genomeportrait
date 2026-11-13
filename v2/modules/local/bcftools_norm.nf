process BCFTOOLS_NORM {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bcftools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.21--h8b25389_0' :
        'biocontainers/bcftools:1.21--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path  fasta

    output:
    tuple val(meta), path("*.norm.vcf.gz"), path("*.norm.vcf.gz.tbi"), emit: vcf
    path "versions.yml",                                               emit: versions

    script:
    def args = task.ext.args ?: '--multiallelics -both'
    """
    bcftools norm $args --fasta-ref $fasta --threads $task.cpus \\
        -Oz -o ${meta.id}.norm.vcf.gz $vcf
    bcftools index --tbi ${meta.id}.norm.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
