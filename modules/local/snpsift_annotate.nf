process SNPSIFT_ANNOTATE {
    tag "$meta.id|$dbtag"
    label 'process_medium'

    // Implemented with bcftools annotate: the SnpSift biocontainer lacks bgzip/tabix, so it
    // could not index ClinVar. bcftools has the full htslib toolchain. Also renames ClinVar's
    // '1'..'MT' contigs to chr-prefixed so they match our GATK hg38 VCF.
    conda "bioconda::bcftools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.21--h8b25389_0' :
        'biocontainers/bcftools:1.21--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path  database
    val   dbtag

    output:
    tuple val(meta), path("*.${dbtag}.vcf.gz"), path("*.${dbtag}.vcf.gz.tbi"), emit: vcf
    path "versions.yml",                                                       emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def fields = dbtag == 'gnomad' ? 'INFO/AF' : (dbtag == 'clinvar' ? 'INFO/CLNSIG,INFO/CLNDN,INFO/CLNREVSTAT' : 'INFO')
    """
    # Map NCBI-style contigs ('1'..'22','X','Y','MT') to chr-prefixed (no-op for chr-prefixed DBs)
    for c in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y; do echo "\$c chr\$c"; done > chr_map.txt
    echo "MT chrM" >> chr_map.txt

    bcftools annotate --rename-chrs chr_map.txt -Oz -o db.${dbtag}.vcf.gz $database
    tabix -p vcf db.${dbtag}.vcf.gz

    bcftools annotate -a db.${dbtag}.vcf.gz -c ${fields} --threads $task.cpus \\
        -Oz -o ${meta.id}.${dbtag}.vcf.gz $vcf
    tabix -p vcf ${meta.id}.${dbtag}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
