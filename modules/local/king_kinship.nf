process KING_KINSHIP {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::king=2.3.2 bioconda::plink2=2.00a5.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/king:2.3.2--h4ac6f70_0' :
        'biocontainers/king:2.3.2--h4ac6f70_0' }"

    input:
    tuple val(meta), path(pgen), path(pvar), path(psam)

    output:
    tuple val(meta), path("*.kin.tsv"), emit: kin
    path "versions.yml",                emit: versions

    script:
    """
    plink2 --pfile ${pgen.baseName} --make-bed --out merged
    king -b merged.bed --kinship --prefix ${meta.id} || true
    # Keep pairs involving the query sample
    if [ -f ${meta.id}.kin0 ]; then
        head -n1 ${meta.id}.kin0 > ${meta.id}.kin.tsv
        grep -w "${meta.id}" ${meta.id}.kin0 >> ${meta.id}.kin.tsv || true
    else
        echo -e "ID1\\tID2\\tKinship\\tnote\\nno_relatives_found" > ${meta.id}.kin.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        king: \$(king 2>&1 | grep -oP 'KING \\K[0-9.]+' | head -n1)
    END_VERSIONS
    """
}
