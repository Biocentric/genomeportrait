process SAMTOOLS_FAIDX {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::samtools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.fai"), emit: fai
    path "versions.yml",            emit: versions

    script:
    // --fai-idx pins the index to the task work dir. Without it samtools writes
    // "<fasta>.fai" beside the resolved source, i.e. INTO the reference storeDir — and since
    // DOWNLOAD_RESOURCE emits a glob ("<id>/*"), the next run then matches BOTH the fasta and
    // the stray .fai, so "$fasta" expands to two paths and faidx reads the .fai as a region.
    """
    samtools faidx --fai-idx ${fasta}.fai $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
