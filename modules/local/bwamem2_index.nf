process BWAMEM2_INDEX {
    tag "$meta.id"
    label 'process_high'
    label 'process_high_memory'
    storeDir "${params.reference_base}/${params.genome}/bwamem2_index"

    conda "bioconda::bwa-mem2=2.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bwa-mem2:2.2.1--he513fc3_0' :
        'biocontainers/bwa-mem2:2.2.1--he513fc3_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    path "bwamem2", emit: index
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p bwamem2
    bwa-mem2 index -p bwamem2/${fasta} $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1 | head -n1)
    END_VERSIONS
    """
}
