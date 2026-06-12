process SNPSIFT_ANNOTATE {
    tag "$meta.id|$dbtag"
    label 'process_medium'

    conda "bioconda::snpsift=5.2 bioconda::htslib=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snpsift:5.2--hdfd78af_0' :
        'biocontainers/snpsift:5.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path  database
    val   dbtag

    output:
    tuple val(meta), path("*.${dbtag}.vcf.gz"), path("*.${dbtag}.vcf.gz.tbi"), emit: vcf
    path "versions.yml",                                                       emit: versions

    script:
    // Pull a sensible set of INFO fields per database
    def info = dbtag == 'gnomad' ? '-info AF,AF_grpmax,nhomalt' : (dbtag == 'clinvar' ? '-info CLNSIG,CLNDN,CLNREVSTAT' : '')
    """
    [ -f "${database}.tbi" ] || tabix -p vcf $database || true
    SnpSift annotate $info $database $vcf | bgzip -c > ${meta.id}.${dbtag}.vcf.gz
    tabix -p vcf ${meta.id}.${dbtag}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpsift: \$(SnpSift 2>&1 | grep -oP 'SnpSift version \\K[0-9.]+' | head -n1)
    END_VERSIONS
    """
}
