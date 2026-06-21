process CRAM_TO_BAM {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(cram), path(crai)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("*.bam"), path("*.bam.bai"), emit: bam
    path "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samtools view -@ $task.cpus -b -T $fasta -o ${meta.id}.bam $cram
    samtools index -@ $task.cpus ${meta.id}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
