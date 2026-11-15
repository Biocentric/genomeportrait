process VERIFYBAMID2 {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::verifybamid2=2.0.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/verifybamid2:2.0.3--hc004090_0' :
        'biocontainers/verifybamid2:2.0.3--hc004090_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("${meta.id}.selfSM"),   emit: selfsm
    tuple val(meta), path("${meta.id}.Ancestry"), emit: ancestry, optional: true
    path "versions.yml",                          emit: versions

    script:
    // DNA from a second person mixed into the library (FREEMIX), estimated jointly with the
    // sample's position on 1000 Genomes PCs, so ancestry is not mistaken for contamination.
    // The GRCh38 marker panel (100k common SNPs) ships inside the bioconda package.
    def panel = params.contamination_panel ?: '/usr/local/share/verifybamid2-2.0.3-0/resource/1000g.phase3.100k.b38.vcf.gz.dat'
    """
    verifybamid2 --SVDPrefix ${panel} --Reference ${fasta} --BamFile ${bam} \\
        --NumThread ${task.cpus} --Output ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        verifybamid2: \$(verifybamid2 --help 2>&1 | grep -o 'Version:[0-9.]*' | cut -d: -f2)
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.id}.selfSM"
    touch "${meta.id}.Ancestry"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: true
    END_VERSIONS
    """
}
