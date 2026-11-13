process GATK4_APPLYBQSR {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::gatk4=4.6.1.0 bioconda::samtools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-d9e7bad0f7fbc8f4458d5c3ab7ffaaf0235b59fb:f857e2d6cc88d35580d01cf39e0959a68b83c1d9-0' :
        'biocontainers/mulled-v2-d9e7bad0f7fbc8f4458d5c3ab7ffaaf0235b59fb:f857e2d6cc88d35580d01cf39e0959a68b83c1d9-0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(table)
    path  fasta
    path  fai
    path  dict

    output:
    tuple val(meta), path("*.recal.cram"), path("*.recal.cram.crai"), emit: cram
    path "versions.yml",                                              emit: versions

    script:
    def avail_mem = (task.memory.mega * 0.8).intValue()
    """
    gatk --java-options "-Xmx${avail_mem}M" ApplyBQSR \\
        --input $bam \\
        --reference $fasta \\
        --bqsr-recal-file $table \\
        --output ${meta.id}.recal.bam
    samtools view -@ $task.cpus -C -T $fasta -o ${meta.id}.recal.cram ${meta.id}.recal.bam
    samtools index ${meta.id}.recal.cram

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -oP 'GATK \\K[0-9.]+')
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
