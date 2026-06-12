process MANTA_GERMLINE {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::manta=1.6.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/manta:1.6.0--h9ee0642_1' :
        'biocontainers/manta:1.6.0--h9ee0642_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("*.diploidSV.vcf.gz"), path("*.diploidSV.vcf.gz.tbi"), emit: sv_vcf
    path "versions.yml",                                                         emit: versions

    script:
    """
    configManta.py --bam $bam --referenceFasta $fasta --runDir manta
    python manta/runWorkflow.py -j $task.cpus
    mv manta/results/variants/diploidSV.vcf.gz     ${meta.id}.diploidSV.vcf.gz
    mv manta/results/variants/diploidSV.vcf.gz.tbi ${meta.id}.diploidSV.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        manta: \$(configManta.py --version)
    END_VERSIONS
    """
}
