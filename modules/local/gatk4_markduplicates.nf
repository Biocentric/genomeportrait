process GATK4_MARKDUPLICATES {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::gatk4=4.6.1.0 bioconda::samtools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-d9e7bad0f7fbc8f4458d5c3ab7ffaaf0235b59fb:f857e2d6cc88d35580d01cf39e0959a68b83c1d9-0' :
        'biocontainers/mulled-v2-d9e7bad0f7fbc8f4458d5c3ab7ffaaf0235b59fb:f857e2d6cc88d35580d01cf39e0959a68b83c1d9-0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("*.markdup.bam"), path("*.markdup.bam.bai"), emit: bam
    tuple val(meta), path("*.metrics"),                                emit: metrics
    path "versions.yml",                                               emit: versions

    script:
    def args = task.ext.args ?: ''
    def avail_mem = (task.memory.mega * 0.8).intValue()
    """
    gatk --java-options "-Xmx${avail_mem}M" MarkDuplicates \\
        --INPUT $bam \\
        --OUTPUT ${meta.id}.markdup.bam \\
        --METRICS_FILE ${meta.id}.markdup.metrics \\
        --REFERENCE_SEQUENCE $fasta \\
        $args
    samtools index ${meta.id}.markdup.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -oP 'GATK \\K[0-9.]+')
    END_VERSIONS
    """
}
