process BLOOD_GROUPS {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? (params.report_sif ?: "${params.reference_base}/containers/genomeportrait-report-${params.report_image_version}.sif") : params.report_docker }"

    input:
    tuple val(meta), path(gvcf), path(tbi), path(rhd)
    path  sites       // assets/blood_groups/sites.pvar
    path  engine      // pgs_engine.py: genotypes at fixed sites from the gVCF, hom-ref included
    path  script

    output:
    tuple val(meta), path("${meta.id}.blood_groups.tsv"), emit: tsv
    path "versions.yml",                                  emit: versions

    script:
    """
    python3 ${engine} target --gvcf ${gvcf} --pvar ${sites} --sample ${meta.id} \\
        --out ${meta.id}.bg.vcf.gz --stats bg_stats.tsv
    python3 ${script} --vcf ${meta.id}.bg.vcf.gz --rhd ${rhd} --sample ${meta.id} --out ${meta.id}.blood_groups.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
