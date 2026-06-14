process PARABRICKS_DEEPVARIANT {
    tag "$meta.id"
    label 'process_high'
    label 'process_gpu'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.6.0-1"

    input:
    tuple val(meta), path(input), path(input_index), path(interval_file)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml",               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                  = task.ext.args ?: ''
    def prefix                = task.ext.prefix ?: "${meta.id}"
    def interval_file_command = interval_file ? interval_file.collect { "--interval-file $it" }.join(' ') : ''
    def num_gpus              = task.accelerator ? "--num-gpus ${task.accelerator.request}" : ''
    """
    pbrun deepvariant \\
        --ref ${fasta} \\
        --in-bam ${input} \\
        ${interval_file_command} \\
        --out-variants ${prefix}.deepvariant.vcf.gz \\
        ${num_gpus} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbrun: \$(echo \$(pbrun version 2>&1) | sed 's/^Please.* //; s/ .*\$//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.deepvariant.vcf.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbrun: 4.6.0-1
    END_VERSIONS
    """
}
