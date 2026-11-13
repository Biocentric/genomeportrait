process TELSEQ {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::telseq=0.0.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/telseq:0.0.2--h06902ac_8' :
        'biocontainers/telseq:0.0.2--h06902ac_8' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.telseq.tsv"), emit: tsv
    path "versions.yml",                   emit: versions

    script:
    """
    telseq -o ${meta.id}.telseq.tsv $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        telseq: \$(telseq --version 2>&1 | grep -oP '[0-9.]+' | head -n1)
    END_VERSIONS
    """
}
