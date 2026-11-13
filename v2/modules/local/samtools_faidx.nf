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
    // NOTE: deliberately plain `samtools faidx <fasta>`.
    // samtools resolves the staged symlink and writes "<fasta>.fai" beside the real file, i.e.
    // into the DOWNLOAD_RESOURCE storeDir. Since that module emits a glob ("<id>/*"), a later
    // run would then match both the FASTA and the stray .fai and expand $fasta to two paths.
    // That is neutralised in prepare_genome.nf, which selects the FASTA out of the glob.
    // Do NOT "fix" it here with --fai-idx: changing this script changes the task hash, which
    // produces a new .fai and cascades a cache miss through markdup/BQSR/DeepVariant/VEP —
    // i.e. a full multi-hour re-run. Only change it when a full re-run is acceptable anyway.
    """
    samtools faidx $fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
