process QUALIMAP_BAMQC {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::qualimap=2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/qualimap:2.3--hdfd78af_0' :
        'biocontainers/qualimap:2.3--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}_qualimap"), emit: results
    path "versions.yml",                          emit: versions

    script:
    def avail_mem = (task.memory.giga).intValue()
    """
    unset DISPLAY
    qualimap --java-mem-size=${avail_mem}G bamqc \\
        -bam $bam \\
        -outdir ${meta.id}_qualimap \\
        -nt $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qualimap: \$(qualimap --help 2>&1 | grep -oP 'QualiMap v\\.\\K[0-9.]+' | head -n1)
    END_VERSIONS
    """
}
