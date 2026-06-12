process PHARMCAT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::pharmcat=2.15.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pharmcat:2.15.5--hdfd78af_0' :
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
    # Preprocess to PharmCAT-ready VCF, then run the full pipeline
    pharmcat_pipeline $vcf \\
        -o pharmcat_out \\
        -bf ${meta.id} \\
        --reference $fasta || true
    cp pharmcat_out/${meta.id}*.report.json ${meta.id}.report.json
    cp pharmcat_out/${meta.id}*.report.html ${meta.id}.report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pharmcat: \$(echo 2.15.5)
    END_VERSIONS
    """
}
