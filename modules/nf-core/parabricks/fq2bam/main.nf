process PARABRICKS_FQ2BAM {
    tag "$meta.id"
    label 'process_high'
    label 'process_gpu'
    stageInMode 'copy'

    // NVIDIA Clara Parabricks — GPU container only (no conda/singularity-from-bioconda).
    // Pulled by Singularity/Apptainer directly from the NVIDIA registry.
    container "${params.parabricks_container}"

    input:
    tuple val(meta), path(reads), path(interval_file)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(index)            // classic BWA index dir (.amb/.ann/.bwt/.pac/.sa)
    tuple val(meta4), path(known_sites)      // optional VCF(s) for the BQSR table

    output:
    tuple val(meta), path("*.bam")                  , emit: bam,               optional:true
    tuple val(meta), path("*.cram")                 , emit: cram,              optional:true
    tuple val(meta), path("*.bai")                  , emit: bai,               optional:true
    tuple val(meta), path("*.table")                , emit: bqsr_table,        optional:true
    tuple val(meta), path("*qc_metrics")            , emit: qc_metrics,        optional:true
    tuple val(meta), path("*.duplicate-metrics.txt"), emit: duplicate_metrics, optional:true
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                  = task.ext.args ?: ''
    def prefix                = task.ext.prefix ?: "${meta.id}"
    // Build --in-fq pairs from the (paired) reads list
    def in_fq_command         = meta.single_end ? "--in-se-fq ${reads}" : "--in-fq ${reads}"
    def known_sites_command   = known_sites ? known_sites.collect { "--knownSites $it" }.join(' ') + " --out-recal-file ${prefix}.table" : ''
    def interval_file_command = interval_file ? interval_file.collect { "--interval-file $it" }.join(' ') : ''
    def num_gpus              = task.accelerator ? "--num-gpus ${task.accelerator.request}" : ''
    """
    # Parabricks' GPU-BWA takes the BWA index prefix from --ref, i.e. it opens
    # <ref>.{amb,ann,bwt,pac,sa} beside the reference. BWA_INDEX emits them in a
    # bwa/ subdirectory under the bare basename ("Homo_sapiens_assembly38.amb",
    # not "Homo_sapiens_assembly38.fasta.amb"), so pbrun would not find them.
    # stageInMode is 'copy', so these are task-local files: move them into place.
    INDEX=\$(find -L . -name "*.amb" | head -n1 | sed 's/\\.amb\$//')
    if [ -n "\${INDEX}" ] && [ "\${INDEX#./}" != "${fasta}" ]; then
        for e in amb ann bwt pac sa; do
            [ -e "${fasta}.\${e}" ] || mv -f "\${INDEX}.\${e}" "${fasta}.\${e}"
        done
    fi

    pbrun fq2bam \\
        --ref ${fasta} \\
        ${in_fq_command} \\
        --read-group-sm ${meta.id} \\
        ${known_sites_command} \\
        ${interval_file_command} \\
        --out-bam ${prefix}.bam \\
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
    touch ${prefix}.bam
    touch ${prefix}.bam.bai
    touch ${prefix}.table
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbrun: 4.6.0-1
    END_VERSIONS
    """
}
