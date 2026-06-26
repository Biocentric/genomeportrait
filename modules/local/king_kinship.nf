process KING_KINSHIP {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::king=2.3.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/king:2.3.2--h4ac6f70_0' :
        'biocontainers/king:2.3.2--h4ac6f70_0' }"

    input:
    tuple val(meta), path(c_bed), path(c_bim), path(c_fam)

    output:
    tuple val(meta), path("*.kin.tsv"), emit: kin
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    king -b ${meta.id}.combined.bed --kinship --prefix ${meta.id} || true
    # Keep pairs involving the query sample
    if [ -f ${meta.id}.kin0 ]; then
        head -n1 ${meta.id}.kin0 > ${meta.id}.kin.tsv
        grep -w "${meta.id}" ${meta.id}.kin0 >> ${meta.id}.kin.tsv || true
    else
        printf "ID1\\tID2\\tKinship\\tnote\\nno_relatives_found\\n" > ${meta.id}.kin.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        king: \$(king 2>&1 | grep -oP 'KING \\K[0-9.]+' | head -n1)
    END_VERSIONS
    """
}
