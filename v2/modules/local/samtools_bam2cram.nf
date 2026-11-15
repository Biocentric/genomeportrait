process SAMTOOLS_BAM2CRAM {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("*.cram"), path("*.cram.crai"), emit: cram
    path "versions.yml",                                  emit: versions

    script:
    // The deliverable alignment is a CRAM (~45 GB at 40x) rather than the duplicate-marked
    // BAM (~120 GB). The output name carries its own suffix so it can never equal the
    // staged input's name: writing through a staged symlink truncates the upstream file.
    def prefix = task.ext.prefix ?: "${meta.id}.md"
    """
    samtools view -@ $task.cpus -C -T $fasta -o ${prefix}.cram $bam
    samtools index -@ $task.cpus ${prefix}.cram

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.id}.cram"
    touch "${meta.id}.cram.crai"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: true
    END_VERSIONS
    """
}
