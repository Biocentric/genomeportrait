process PGS_REPORT {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.12 conda-forge::numpy=2.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? params.report_sif : params.report_docker }"

    input:
    tuple val(meta), path(target_scores), path(target_pcs)
    path  ref_scores
    path  ref_pcs
    tuple path(ref_pgen), path(ref_pvar), path(ref_psam)
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
    tuple val(meta), path("${meta.id}.pgs_interpreted.tsv"),  emit: tsv
    tuple val(meta), path("${meta.id}.pgs_interpreted.json"), emit: json
    tuple val(meta), path("${meta.id}.pgs_overview.html"),    emit: html
    tuple val(meta), path("${meta.id}.pgs_fragment.html"),    emit: fragment
    path "versions.yml",                                      emit: versions

    script:
    """
    python3 ${engine} normalise \\
        --ref-scores ${ref_scores} --target-scores ${target_scores} \\
        --ref-pcs ${ref_pcs} --target-pcs ${target_pcs} \\
        --psam ${ref_psam} --match-stats ${match_stats} \\
        --min-overlap ${params.pgs_min_overlap} \\
        --out ${meta.id}.pgs_scores.tsv --pop-out ${meta.id}.pgs_population.tsv

    python3 ${interpret} --scores ${meta.id}.pgs_scores.tsv --panel ${panel} --concepts ${concepts} \\
        --catalog ${catalog} --r2-exclusions ${r2_exclusions} --sample ${meta.id} \\
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
