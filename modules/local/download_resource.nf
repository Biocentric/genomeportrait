process DOWNLOAD_RESOURCE {
    tag "$meta.id"
    label 'process_single'
    storeDir storedir   // <-- cache: skip if already present, never re-download

    conda "conda-forge::wget=1.21.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wget:1.21.4' :
        'biocontainers/wget:1.21.4' }"

    input:
    tuple val(meta), val(url)
    val storedir

    output:
    tuple val(meta), path("${meta.id}/*"), emit: file
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def fname = url.tokenize('/').last()
    """
    # Download into a per-resource subdir so the original filename (and its extension,
    # which downstream tools rely on) is preserved without fragile shell renaming.
    mkdir -p ${meta.id}
    wget --quiet --continue --tries=3 -O "${meta.id}/${fname}" "${url}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | head -n1 | sed 's/^GNU Wget //; s/ .*//')
    END_VERSIONS
    """
}
