process YLEAF {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::yleaf=3.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/yleaf:3.2.1--pyh1286868_0' :
        'biocontainers/yleaf:3.2.1--pyh1286868_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("*.y_haplogroup.txt"), emit: hg
    path "versions.yml",                         emit: versions

    script:
    """
    Yleaf -bam $bam -rg hg38 -o ${meta.id}_yleaf -t $task.cpus || true
    # Collect the predicted haplogroup
    if [ -f ${meta.id}_yleaf/hg_prediction.hg ]; then
        cp ${meta.id}_yleaf/hg_prediction.hg ${meta.id}.y_haplogroup.txt
    else
        find ${meta.id}_yleaf -name '*.hg' -exec cp {} ${meta.id}.y_haplogroup.txt \\; 2>/dev/null || \\
        echo -e "Sample\\tHaplogroup\\nNA\\tNA(no_Y_reads)" > ${meta.id}.y_haplogroup.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yleaf: \$(Yleaf --version 2>&1 | grep -oP '[0-9.]+' | head -n1)
    END_VERSIONS
    """
}
