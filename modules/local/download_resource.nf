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
    tuple val(meta), path("${meta.id}.*"), emit: file
    path "versions.yml",                   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def fname = url.tokenize('/').last()
    """
    # Download with resume; index bgzipped VCFs that ship without a .tbi
    wget --quiet --continue --tries=3 -O "${fname}" "${url}"

    # Keep meta.id prefix so downstream globbing is deterministic
    if [ "${fname}" != "${meta.id}.${fname#*.}" ]; then
        ln -s "${fname}" "${meta.id}.${fname#*.}" 2>/dev/null || cp "${fname}" "${meta.id}.${fname#*.}"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | head -n1 | sed 's/^GNU Wget //; s/ .*//')
    END_VERSIONS
    """
}
