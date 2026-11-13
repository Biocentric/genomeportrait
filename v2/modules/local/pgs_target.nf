process PGS_TARGET {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.12 conda-forge::numpy=2.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? params.report_sif : params.report_docker }"

    input:
    tuple val(meta), path(gvcf), path(tbi)
    tuple path(ref_pgen), path(ref_pvar), path(ref_psam)
    path  engine

    output:
    tuple val(meta), path("${meta.id}.pgs_target.vcf.gz"),    emit: vcf
    tuple val(meta), path("${meta.id}.pgs_target_stats.tsv"), emit: stats
    path "versions.yml",                                      emit: versions

    script:
    // A score is a weighted sum over FIXED positions, so every scored position needs a
    // genotype. The variant-only VCF simply omits hom-ref sites; the gVCF asserts them in
    // reference blocks. Rebuilding the sample in the reference panel's allele space (rather
    // than expanding the gVCF, which keeps its <*> allele) is what lets alleles match.
    """
    python3 ${engine} target --gvcf ${gvcf} --pvar ${ref_pvar} --sample ${meta.id} \\
        --min-gq ${params.pgs_min_gq} \\
        --out ${meta.id}.pgs_target.vcf.gz --stats ${meta.id}.pgs_target_stats.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
