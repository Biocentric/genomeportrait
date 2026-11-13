process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::minimap2=2.28 bioconda::samtools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' :
        'biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' }"

    input:
    tuple val(meta), path(reads)
    path  fasta

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml",            emit: versions

    script:
    def args = task.ext.args ?: ''
    // ONT vs PacBio HiFi presets
    def preset = meta.platform == 'pacbio' ? 'map-hifi' : 'map-ont'
    """
    minimap2 -ax $preset $args -t $task.cpus $fasta ${reads} \\
        | samtools sort -@ $task.cpus -o ${meta.id}.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
    END_VERSIONS
    """
}
