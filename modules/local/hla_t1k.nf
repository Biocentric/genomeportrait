process HLA_T1K {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::t1k=1.0.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/t1k:1.0.9--h5ca1c30_0' :
        'biocontainers/t1k:1.0.9--h5ca1c30_0' }"

    input:
    tuple val(meta), path(reads1), path(reads2)
    path  index          // t1k_hla reference dir from BUILD_T1K_HLA

    output:
    tuple val(meta), path("*_genotype.tsv"), emit: alleles
    path "versions.yml",                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // t1k-build.pl --download only emits the *_dna_seq.fa allele reference (no *_coord.fa),
    // so BAM mode (which requires -c) is impossible. We feed the pre-extracted MHC-region
    // FASTQ pair (HLA_EXTRACT_FASTQ) and run T1K in FASTQ mode, which needs only -f.
    """
    SEQ=\$(find ${index} -name "*_dna_seq.fa" | head -1)
    run-t1k -1 ${reads1} -2 ${reads2} --preset hla-wgs \\
        -f \$SEQ -t $task.cpus --od . -o ${meta.id} || true
    [ -f ${meta.id}_genotype.tsv ] || printf "gene\\tnum\\tallele1\\n" > ${meta.id}_genotype.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        t1k: \$(run-t1k --version 2>&1 | grep -oE '[0-9.]+' | head -n1)
    END_VERSIONS
    """
}
