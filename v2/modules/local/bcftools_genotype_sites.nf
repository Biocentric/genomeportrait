process BCFTOOLS_GENOTYPE_SITES {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bcftools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.21--h8b25389_0' :
        'biocontainers/bcftools:1.21--h8b25389_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fai
    path  sites          // CHROM POS REF ALT (chr-prefixed)

    output:
    tuple val(meta), path("${meta.id}.anc.vcf.gz"), path("${meta.id}.anc.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Genotype the sample at the fixed ancestry SNP positions directly from the alignment,
    // WITHOUT -v, so homozygous-reference sites are emitted (the DeepVariant VCF is
    // variant-biased and gives nonsense ancestry).
    """
    cut -f1,2 $sites > regions.tsv
    bcftools mpileup -R regions.tsv -f $fasta -a FORMAT/DP -Ou $bam 2>/dev/null \\
        | bcftools call -m -Oz -o ${meta.id}.anc.vcf.gz
    bcftools index -t ${meta.id}.anc.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
