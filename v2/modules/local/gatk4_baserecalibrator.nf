process GATK4_BASERECALIBRATOR {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::gatk4=4.6.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.6.1.0--py310hdfd78af_0' :
        'biocontainers/gatk4:4.6.1.0--py310hdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fai
    path  dict
    path  dbsnp
    path  known_indels
    path  mills

    output:
    tuple val(meta), path("*.recal.table"), emit: table
    path "versions.yml",                    emit: versions

    script:
    def avail_mem = (task.memory.mega * 0.8).intValue()
    def sites = "--known-sites ${dbsnp} --known-sites ${known_indels} --known-sites ${mills}"
    """
    # Ensure known-site VCFs are tabix-indexed
    for v in ${dbsnp} ${known_indels} ${mills}; do
        [ -f "\${v}.tbi" ] || [ -f "\${v}.idx" ] || gatk IndexFeatureFile -I "\${v}"
    done

    gatk --java-options "-Xmx${avail_mem}M" BaseRecalibrator \\
        --input $bam \\
        --reference $fasta \\
        $sites \\
        --output ${meta.id}.recal.table

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -oP 'GATK \\K[0-9.]+')
    END_VERSIONS
    """
}
