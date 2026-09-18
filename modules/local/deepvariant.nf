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
    // make_examples writes tens of GB of TFRecords from every shard at once. Left in the
    // task directory those land on whatever filesystem the work dir uses - on a FUSE
    // mount (ntfs-3g) that single userspace process becomes the serialisation point for
    // every shard. Point this at a local kernel filesystem (e.g. an NVMe scratch bound
    // into the container) when the work dir is not one.
    def interdir = params.dv_intermediate_dir ? "${params.dv_intermediate_dir}/dv_${meta.id}" : 'tmp'
    """
    mkdir -p ${interdir}
    trap 'rm -rf ${interdir}' EXIT

    run_deepvariant \\
        $args \\
        --ref=$fasta \\
        --reads=$bam \\
        --output_vcf=${meta.id}.deepvariant.vcf.gz \\
        --output_gvcf=${meta.id}.g.vcf.gz \\
        --num_shards=$task.cpus \\
        --intermediate_results_dir=${interdir}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        deepvariant: \$(run_deepvariant --version 2>&1 | grep -oP '[0-9]+\\.[0-9]+\\.[0-9]+' | head -n1)
    END_VERSIONS
    """
}
