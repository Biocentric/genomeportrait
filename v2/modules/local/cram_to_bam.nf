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
    tuple val(meta), path("*.tobam.bam"), path("*.tobam.bam.bai"), emit: bam
    path "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // The output must not collide with the staged input name. When the analysis alignment
    // is already a BAM (Parabricks path, or --skip_bqsr) the input is staged as a symlink
    // called <meta.id>.bam, and writing an output of the same name truncates the upstream
    // file through that symlink. The distinct suffix also keeps the output glob from
    // matching the staged input.
    def prefix = task.ext.prefix ?: "${meta.id}.tobam"
    """
    samtools view -@ $task.cpus -b -T $fasta -o ${prefix}.bam $cram
    samtools index -@ $task.cpus ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
