process TRAIT_GENOTYPE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bcftools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.21--h8b25389_0' :
        'biocontainers/bcftools:1.21--h8b25389_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fai
    path  panel

    output:
    tuple val(meta), path("*.trait_geno.tsv"), emit: geno
    path "versions.yml",                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Targeted genotyping at the trait-panel positions, so homozygous-reference sites are
    # captured (DeepVariant emits variant sites only). Build BED from the panel coords.
    awk -F'\\t' '\$1 !~ /^#/ && \$10 ~ /^chr/ && \$11 ~ /^[0-9]+\$/ {print \$10"\\t"(\$11-1)"\\t"\$11}' $panel \\
        | sort -k1,1 -k2,2n -u > trait.bed

    bcftools mpileup -R trait.bed -f $fasta -a FORMAT/DP $bam -Ou 2>/dev/null \\
        | bcftools call -m -Oz -o calls.vcf.gz
    bcftools index -t calls.vcf.gz

    # One row per covered panel position: CHROM POS REF ALT GT DP
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT[\\t%GT\\t%DP]\\n' calls.vcf.gz > ${meta.id}.trait_geno.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
