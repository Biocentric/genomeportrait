process MUTSERVE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::mutserve=2.0.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mutserve:2.0.3--hdfd78af_0' :
        'biocontainers/mutserve:2.0.3--hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("*.heteroplasmy.txt"), emit: txt
    path "versions.yml",                          emit: versions

    script:
    """
    # Detects low-level mtDNA heteroplasmies (variant level >= 1%)
    mutserve call \\
        --reference $fasta \\
        --output ${meta.id}.heteroplasmy.txt \\
        --level 0.01 \\
        $bam || echo -e "Pos\\tRef\\tVariant\\tVariantLevel" > ${meta.id}.heteroplasmy.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mutserve: \$(mutserve version 2>&1 | grep -oP '[0-9.]+' | head -n1)
    END_VERSIONS
    """
}
