process VCF_TO_TABLE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bcftools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.21--h8b25389_0' :
        'biocontainers/bcftools:1.21--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.annotation.tsv.gz"), emit: tsv
    path "versions.yml",                          emit: versions

    script:
    """
    # Flatten VEP CSQ + key INFO fields into one row per consequence
    bcftools +split-vep -d -f \\
      '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%SYMBOL\\t%Consequence\\t%IMPACT\\t%gnomADg_AF\\t%CLNSIG\\t%Existing_variation\\t[%GT]\\n' \\
      -A tab $vcf 2>/dev/null \\
      | sed '1i CHROM\\tPOS\\tID\\tREF\\tALT\\tGENE\\tCONSEQUENCE\\tIMPACT\\tGNOMAD_AF\\tCLINSIG\\tDBSNP\\tGT' \\
      | bgzip -c > ${meta.id}.annotation.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
