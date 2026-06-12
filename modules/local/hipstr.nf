process HIPSTR {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::hipstr=0.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hipstr:0.7--he4a0461_1' :
        'biocontainers/hipstr:0.7--he4a0461_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fai
    path  str_bed

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
