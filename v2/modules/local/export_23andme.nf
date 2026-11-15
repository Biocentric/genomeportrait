process EXPORT_23ANDME {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? (params.report_sif ?: "${params.reference_base}/containers/genomeportrait-report-${params.report_image_version}.sif") : params.report_docker }"

    input:
    tuple val(meta), path(gvcf), path(tbi), val(sex)
    path  sites
    path  fasta
    path  fai
    path  script
    path  engine      // pgs_engine.py: tabix/BGZF reading shared with the PGS step

    output:
    tuple val(meta), path("${meta.id}.23andme.txt.gz"), emit: txt
    path "versions.yml",                                emit: versions

    script:
    """
    python3 ${script} --gvcf ${gvcf} --sites ${sites} --fasta ${fasta} --sex ${sex} \\
        --sample ${meta.id} --threads ${task.cpus} --out ${meta.id}.23andme.txt
    gzip -f ${meta.id}.23andme.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
