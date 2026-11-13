process ENSEMBLVEP_VEP {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::ensembl-vep=113.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ensembl-vep:113.0--pl5321h2a3209d_0' :
        'biocontainers/ensembl-vep:113.0--pl5321h2a3209d_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path  fasta
    path  vep_cache
    val   species
    val   cache_version

    output:
    tuple val(meta), path("*.vep.vcf.gz"), path("*.vep.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("*.vep.summary.html"),                     emit: report
    path "versions.yml",                                             emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    vep \\
        --input_file $vcf \\
        --output_file ${meta.id}.vep.vcf.gz \\
        --stats_file ${meta.id}.vep.summary.html \\
        --species $species \\
        --cache_version $cache_version \\
        --dir_cache $vep_cache \\
        --fasta $fasta \\
        --fork $task.cpus \\
        --vcf --compress_output bgzip --force_overwrite \\
        $args
    tabix -p vcf ${meta.id}.vep.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensembl-vep: \$(vep --help 2>&1 | grep -oP 'ensembl-vep\\s+:\\s+\\K[0-9.]+')
    END_VERSIONS
    """
}
