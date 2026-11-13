process GATK4_HAPLOTYPECALLER {
    tag "$meta.id"
    label 'process_high'

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

    output:
    tuple val(meta), path("*.haplotypecaller.vcf.gz"), path("*.haplotypecaller.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("*.g.vcf.gz"), path("*.g.vcf.gz.tbi"),                              emit: gvcf
    path "versions.yml",                                                                      emit: versions

    script:
    def avail_mem = (task.memory.mega * 0.8).intValue()
    """
    gatk --java-options "-Xmx${avail_mem}M" HaplotypeCaller \\
        --input $bam --reference $fasta --dbsnp $dbsnp \\
        --output ${meta.id}.g.vcf.gz -ERC GVCF \\
        --native-pair-hmm-threads $task.cpus
    gatk --java-options "-Xmx${avail_mem}M" GenotypeGVCFs \\
        --reference $fasta --variant ${meta.id}.g.vcf.gz \\
        --output ${meta.id}.haplotypecaller.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(gatk --version 2>&1 | grep -oP 'GATK \\K[0-9.]+')
    END_VERSIONS
    """
}
