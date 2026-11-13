process DOWNLOAD_VEP_CACHE {
    tag "vep_cache"
    label 'process_single'
    label 'process_long'
    storeDir storedir

    conda "conda-forge::wget=1.21.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wget:1.21.4' :
        'biocontainers/wget:1.21.4' }"

    input:
    val url
    val storedir

    output:
    path "vep_cache",    emit: cache
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p vep_cache
    wget --quiet --continue --tries=3 -O vep_cache.tar.gz "${url}"
    tar -xzf vep_cache.tar.gz -C vep_cache
    rm vep_cache.tar.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | head -n1 | sed 's/^GNU Wget //; s/ .*//')
    END_VERSIONS
    """
}
