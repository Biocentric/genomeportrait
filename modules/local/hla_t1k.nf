process HLA_T1K {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::t1k=1.0.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/t1k:1.0.9--h5ca1c30_0' :
        'biocontainers/t1k:1.0.9--h5ca1c30_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fai
    path  index          // t1k_hla reference dir from BUILD_T1K_HLA

    output:
    tuple val(meta), path("*_genotype.tsv"), emit: alleles
    path "versions.yml",                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    SEQ=\$(find ${index} -name "*_dna_seq.fa"   | head -1)
    COORD=\$(find ${index} -name "*_dna_coord.fa" | head -1)
    run-t1k -b $bam --preset hla-wgs \\
        -f \$SEQ -c \$COORD \\
        -t $task.cpus --od . -o ${meta.id} || true
    [ -f ${meta.id}_genotype.tsv ] || printf "gene\\tnum\\tallele1\\n" > ${meta.id}_genotype.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        t1k: \$(run-t1k --version 2>&1 | grep -oE '[0-9.]+' | head -n1)
    END_VERSIONS
    """
}
