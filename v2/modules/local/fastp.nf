process FASTP {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::fastp=0.23.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.4--h5f740d0_0' :
        'biocontainers/fastp:0.23.4--h5f740d0_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.trim.fastq.gz"), emit: reads
    tuple val(meta), path("*.json"),          emit: json
    tuple val(meta), path("*.html"),          emit: html
    path "versions.yml",                      emit: versions

    script:
    def args = task.ext.args ?: ''
    def paired = reads instanceof List && reads.size() == 2
    if (paired) {
        """
        fastp -i ${reads[0]} -I ${reads[1]} \\
            -o ${meta.id}_1.trim.fastq.gz -O ${meta.id}_2.trim.fastq.gz \\
            --json ${meta.id}.fastp.json --html ${meta.id}.fastp.html \\
            --thread $task.cpus $args 2> ${meta.id}.fastp.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed 's/fastp //')
        END_VERSIONS
        """
    } else {
        """
        fastp -i ${reads} -o ${meta.id}.trim.fastq.gz \\
            --json ${meta.id}.fastp.json --html ${meta.id}.fastp.html \\
            --thread $task.cpus $args 2> ${meta.id}.fastp.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed 's/fastp //')
        END_VERSIONS
        """
    }
}
