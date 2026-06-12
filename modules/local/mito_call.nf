process MITO_CALL {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bcftools=1.21 bioconda::samtools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.21--h8b25389_0' :
        'biocontainers/bcftools:1.21--h8b25389_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("*.chrM.vcf"), emit: vcf
    path "versions.yml",                 emit: versions

    script:
    // Detect MT contig name (chrM vs MT) then call a haploid consensus for HaploGrep
    """
    MT=\$(grep -E '^(chrM|MT)\\b' ${fai} | cut -f1 | head -n1)
    bcftools mpileup -f $fasta -r \${MT} $bam \\
        | bcftools call --ploidy 1 -m -Ov \\
        | bcftools norm -f $fasta -Ov -o ${meta.id}.chrM.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
