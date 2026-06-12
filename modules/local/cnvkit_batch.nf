process CNVKIT_BATCH {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::cnvkit=0.9.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvkit:0.9.11--pyhdfd78af_0' :
        'biocontainers/cnvkit:0.9.11--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta

    output:
    tuple val(meta), path("*.cnr"), emit: cnr
    tuple val(meta), path("*.cns"), emit: cns
    path "versions.yml",            emit: versions

    script:
    def args = task.ext.args ?: '--method wgs'
    """
    cnvkit.py batch $bam -n \\
        --fasta $fasta $args \\
        --output-dir cnvkit -p $task.cpus
    cp cnvkit/*.cnr ${meta.id}.cnr
    cp cnvkit/*.cns ${meta.id}.cns

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvkit: \$(cnvkit.py version | sed 's/cnvkit v//')
    END_VERSIONS
    """
}
