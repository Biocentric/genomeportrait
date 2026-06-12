process DEEPVARIANT_RUNDEEPVARIANT {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://google/deepvariant:1.6.1' :
        'google/deepvariant:1.6.1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("*.deepvariant.vcf.gz"),  path("*.deepvariant.vcf.gz.tbi"),  emit: vcf
    tuple val(meta), path("*.g.vcf.gz"), path("*.g.vcf.gz.tbi"),                       emit: gvcf
    path "versions.yml",                                                               emit: versions

    script:
    def args = task.ext.args ?: "--model_type=${params.dv_model_type}"
    """
    run_deepvariant \\
        $args \\
        --ref=$fasta \\
        --reads=$bam \\
        --output_vcf=${meta.id}.deepvariant.vcf.gz \\
        --output_gvcf=${meta.id}.g.vcf.gz \\
        --num_shards=$task.cpus \\
        --intermediate_results_dir=tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant: \$(run_deepvariant --version 2>&1 | grep -oP '[0-9]+\\.[0-9]+\\.[0-9]+' | head -n1)
    END_VERSIONS
    """
}
