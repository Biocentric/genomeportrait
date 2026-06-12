process BWAMEM2_MEM {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' :
        'biocontainers/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' }"

    input:
    tuple val(meta), path(reads)
    path  index
    path  fasta
    val   sort_bam

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml",            emit: versions

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def sort  = sort_bam ? "samtools sort -@ ${task.cpus} -o ${meta.id}.bam -" : "samtools view -@ ${task.cpus} -bhS -o ${meta.id}.bam -"
    """
    INDEX=\$(find -L ./ -name "*.amb" | sed 's/\\.amb\$//')
    bwa-mem2 mem $args -t $task.cpus \$INDEX ${reads} | $sort

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1 | head -n1)
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
