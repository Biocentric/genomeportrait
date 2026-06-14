process CNVPYTOR {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::cnvpytor=1.3.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvpytor:1.3.2--pyhdfd78af_0' :
        'biocontainers/cnvpytor:1.3.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("*.cnvpytor.tsv"), emit: calls
    path "versions.yml",                     emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Read-depth CNV (CNVnator successor) — efficient single-sample germline WGS.
    def bin = task.ext.args ?: '100000'
    """
    cnvpytor -root ${meta.id}.pytor -rd ${bam}
    cnvpytor -root ${meta.id}.pytor -his ${bin}
    cnvpytor -root ${meta.id}.pytor -partition ${bin}
    cnvpytor -root ${meta.id}.pytor -call ${bin} > ${meta.id}.cnvpytor.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvpytor: \$(cnvpytor --version 2>&1 | grep -oP '[0-9.]+' | head -n1)
    END_VERSIONS
    """
}
