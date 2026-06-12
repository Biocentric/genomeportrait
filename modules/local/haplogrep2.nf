process HAPLOGREP2 {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::haplogrep=2.4.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/haplogrep:2.4.0--hdfd78af_0' :
        'biocontainers/haplogrep:2.4.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(mt_vcf)

    output:
    tuple val(meta), path("*.mtdna_haplogroup.txt"), emit: result
    path "versions.yml",                              emit: versions

    script:
    """
    haplogrep classify \\
        --in $mt_vcf --format vcf \\
        --out ${meta.id}.mtdna_haplogroup.txt \\
        --extend-report

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        haplogrep: \$(haplogrep --version 2>&1 | grep -oP '[0-9]+\\.[0-9]+\\.[0-9]+' | head -n1)
    END_VERSIONS
    """
}
