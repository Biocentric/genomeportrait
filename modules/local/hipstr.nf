process HIPSTR {
    tag "$meta.id"
    label 'process_medium'

    // HipSTR is not on bioconda — use the image built at runtime by BUILD_HIPSTR_CONTAINER.
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? params.hipstr_sif : params.hipstr_docker }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fai
    path  str_bed
    val   ready          // gate: ensures BUILD_HIPSTR_CONTAINER finished before this runs

    output:
    tuple val(meta), path("*.hipstr.vcf.gz"), emit: vcf
    path "versions.yml",                      emit: versions

    script:
    """
    REGIONS=${str_bed}
    if [[ ${str_bed} == *.gz ]]; then gunzip -kf ${str_bed}; REGIONS=${str_bed.baseName}; fi
    HipSTR \\
        --bams $bam \\
        --fasta $fasta \\
        --regions \${REGIONS} \\
        --str-vcf ${meta.id}.hipstr.vcf.gz \\
        --output-filters

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hipstr: \$(HipSTR --version 2>&1 | grep -oP '[0-9.]+' | head -n1)
    END_VERSIONS
    """
}
