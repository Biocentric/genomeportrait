process BUILD_T1K_KIR {
    tag "t1k_kir"
    label 'process_medium'
    label 'process_long'
    storeDir "${params.reference_base}/t1k_kir"

    conda "bioconda::t1k=1.0.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/t1k:1.0.9--h5ca1c30_0' :
        'biocontainers/t1k:1.0.9--h5ca1c30_0' }"

    output:
    path "t1k_kir",      emit: index
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p t1k_kir
    # Builds <prefix>_dna_seq.fa + <prefix>_dna_coord.fa from the IPD-KIR database.
    t1k-build.pl -o t1k_kir --download IPD-KIR

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        t1k: \$(run-t1k --version 2>&1 | grep -oE '[0-9.]+' | head -n1)
    END_VERSIONS
    """
}
