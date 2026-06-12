process HLA_T1K {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::t1k=1.0.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/t1k:1.0.5--hdcf5f25_1' :
        'biocontainers/t1k:1.0.5--hdcf5f25_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("*_genotype.tsv"), emit: alleles
    path "versions.yml",                     emit: versions

    script:
    """
    # T1K ships the IPD-IMGT/HLA reference; build coordinate files on first use
    HLA_REF=\$(find / -name "hla.fa" -path "*t1k*" 2>/dev/null | head -n1)
    HLA_COORD=\$(find / -name "hla_coord.fa" -path "*t1k*" 2>/dev/null | head -n1)
    run-t1k \\
        -b $bam \\
        --preset hla-wgs \\
        -f \${HLA_REF} -c \${HLA_COORD} \\
        -t $task.cpus \\
        -o ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        t1k: \$(run-t1k --version 2>&1 | grep -oP '[0-9.]+' | head -n1)
    END_VERSIONS
    """
}
