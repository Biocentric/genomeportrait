process BWA_INDEX {
    tag "$meta.id"
    label 'process_high'
    label 'process_high_memory'
    storeDir "${params.reference_base}/${params.genome}/bwa_index"

    conda "bioconda::bwa=0.7.18"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwa:0.7.18--he4a0461_1' :
        'biocontainers/bwa:0.7.18--he4a0461_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("bwa"), emit: index
    path "versions.yml",          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p bwa
    bwa index -p bwa/${fasta.baseName} ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact.*\$//')
    END_VERSIONS
    """
}
