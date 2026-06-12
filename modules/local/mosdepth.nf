process MOSDEPTH {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::mosdepth=0.3.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mosdepth:0.3.10--h4e814b3_1' :
        'biocontainers/mosdepth:0.3.10--h4e814b3_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta

    output:
    tuple val(meta), path("*.mosdepth.summary.txt"),    emit: summary
    tuple val(meta), path("*.mosdepth.global.dist.txt"),emit: global_dist
    tuple val(meta), path("*.regions.bed.gz"),          emit: regions, optional: true
    path "versions.yml",                                emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    mosdepth --threads $task.cpus --fasta $fasta $args ${meta.id} $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mosdepth: \$(mosdepth --version | sed 's/mosdepth //')
    END_VERSIONS
    """
}
