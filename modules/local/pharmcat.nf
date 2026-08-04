process PHARMCAT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::pharmcat=2.15.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://pgkb/pharmcat:2.15.5' :
        'pgkb/pharmcat:2.15.5' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("*.report.json"), emit: json
    tuple val(meta), path("*.report.html"), emit: html
    path "versions.yml",                    emit: versions

    script:
    """
    # Preprocess to PharmCAT-ready VCF, then run the full pipeline.
    # (pharmcat_pipeline uses its bundled reference positions — no --reference flag.)
    # -reporterJson is REQUIRED: without it PharmCAT writes only the HTML report, the
    # .report.json copy below falls back to "{}", and the downstream reporter then dies with
    # "No columns to parse from file".
    pharmcat_pipeline $vcf -o pharmcat_out -bf ${meta.id} -reporterJson || true
    # Fall back to phenotype.json (always written) if the reporter JSON is missing.
    if [ -s pharmcat_out/${meta.id}.report.json ]; then
        cp pharmcat_out/${meta.id}.report.json ${meta.id}.report.json
    elif [ -s pharmcat_out/${meta.id}.phenotype.json ]; then
        echo "no reporter JSON; falling back to phenotype.json" >&2
        cp pharmcat_out/${meta.id}.phenotype.json ${meta.id}.report.json
    else
        cp pharmcat_out/${meta.id}*.report.json ${meta.id}.report.json 2>/dev/null || echo '{}' > ${meta.id}.report.json
    fi
    cp pharmcat_out/${meta.id}*.report.html ${meta.id}.report.html 2>/dev/null || echo '<html><body>PharmCAT produced no report</body></html>' > ${meta.id}.report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pharmcat: \$(echo 2.15.5)
    END_VERSIONS
    """
}
