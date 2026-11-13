process PGS_REPORT {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.12 conda-forge::numpy=2.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? (params.report_sif ?: "${params.reference_base}/containers/genomeportrait-report-${params.report_image_version}.sif") : params.report_docker }"

    input:
    tuple val(meta), path(target_vcf), path(target_pcs)
    path  ref_scores
    path  ref_pcs
    tuple path(ref_pgen), path(ref_pvar), path(ref_psam)
    path  ref_afreq
    path  pop_afreq
    path  scores
    path  match_stats
    path  panel
    path  concepts
    path  r2_exclusions
    path  catalog
    path  engine
    path  interpret
    path  overview

    output:
    tuple val(meta), path("${meta.id}.pgs_scores.tsv"),       emit: scores
    tuple val(meta), path("${meta.id}.pgs_population.tsv"),   emit: population
    tuple val(meta), path("${meta.id}.pgs_imputation.tsv"),   emit: imputation
    tuple val(meta), path("${meta.id}.pgs_interpreted.tsv"),  emit: tsv
    tuple val(meta), path("${meta.id}.pgs_interpreted.json"), emit: json
    tuple val(meta), path("${meta.id}.pgs_overview.html"),    emit: html
    tuple val(meta), path("${meta.id}.pgs_fragment.html"),    emit: fragment
    path "versions.yml",                                      emit: versions

    script:
    """
    # Most similar reference population: its allele frequencies fill missing genotypes.
    python3 ${engine} population --ref-pcs ${ref_pcs} --target-pcs ${target_pcs} --psam ${ref_psam} \\
        --pop-col ${params.pgs_population_col} --out population_for_imputation.tsv

    # One pass over the sample's genotypes; every score is a dot product over its .npz rows.
    python3 ${engine} score --target ${target_vcf} --pvar ${ref_pvar} --scores ${scores} \\
        --impute-afreq ${ref_afreq} --pop-afreq-dir ${pop_afreq} \\
        --population population_for_imputation.tsv --impute ${params.pgs_impute} \\
        --imputation-out ${meta.id}.pgs_imputation.tsv --outdir target_scores

    python3 ${engine} normalise \\
        --ref-scores ${ref_scores} --target-scores target_scores \\
        --ref-pcs ${ref_pcs} --target-pcs ${target_pcs} \\
        --psam ${ref_psam} --match-stats ${match_stats} \\
        --min-overlap ${params.pgs_min_overlap} \\
        --out ${meta.id}.pgs_scores.tsv --pop-out ${meta.id}.pgs_population.tsv

    python3 ${interpret} --scores ${meta.id}.pgs_scores.tsv --panel ${panel} --concepts ${concepts} \\
        --catalog ${catalog} --r2-exclusions ${r2_exclusions} --sample ${meta.id} \\
        --imputation ${meta.id}.pgs_imputation.tsv \\
        --percentile ${params.pgs_percentile} \\
        --out-tsv ${meta.id}.pgs_interpreted.tsv --out-json ${meta.id}.pgs_interpreted.json

    python3 ${overview} ${meta.id}.pgs_interpreted.json \\
        --html ${meta.id}.pgs_overview.html --fragment ${meta.id}.pgs_fragment.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
