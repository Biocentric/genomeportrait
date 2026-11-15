process BWAMEM2_MARKDUP_CRAM {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::bwa-mem2=2.2.1 bioconda::samtools=1.16.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' :
        'biocontainers/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2cdf6bf1e92acbeb9b2834b1c58754167173a410-0' }"

    input:
    tuple val(meta), path(reads)
    path  index
    path  alt        // bwa .alt file for the reference ([] = none): ALT-aware mapping
    path  fasta
    path  fai

    output:
    tuple val(meta), path("${meta.id}.md.cram"), path("${meta.id}.md.cram.crai"), emit: cram
    tuple val(meta), path("${meta.id}.markdup.stats"),                           emit: metrics
    path "versions.yml",                                                         emit: versions

    script:
    // Reads to a duplicate-marked CRAM in one stream:
    //   bwa-mem2 -> fixmate -m (mate scores, needed by markdup) -> coordinate sort -> markdup -> CRAM
    // The only intermediates are sort's temporary files. The separate path (sorted BAM, then
    // GATK MarkDuplicates BAM, then CRAM) writes two whole-genome BAMs, ~240 GB at 40x.
    // Uncompressed BAM (-u) between the steps avoids compressing data that is read again at once.
    def args   = task.ext.args ?: ''
    def sort_threads = Math.max(1, (task.cpus / 4) as int)
    def mem_per_thread = task.memory ? "${Math.max(768, ((task.memory.toMega() * 0.25) / sort_threads) as int)}M" : '768M'
    """
    INDEX=\$(find -L ./ -name "*.amb" | sed 's/\\.amb\$//')
    # ALT-aware mapping: bwa-mem2 reads <index prefix>.alt. The cached index directory is not
    # written to; its files are linked into the task directory next to the .alt file instead.
    if [ -s "${alt}" ]; then
        mkdir -p altidx
        for f in \${INDEX}.*; do ln -sf "\$(readlink -f \$f)" altidx/; done
        ln -sf "\$(readlink -f ${alt})" "altidx/\$(basename \$INDEX).alt"
        INDEX=altidx/\$(basename \$INDEX)
    fi
    mkdir -p tmp
    bwa-mem2 mem $args -t $task.cpus \$INDEX ${reads} \\
        | samtools fixmate -m -u - - \\
        | samtools sort -u -@ ${sort_threads} -m ${mem_per_thread} -T tmp/sort - \\
        | samtools markdup -@ ${sort_threads} -T tmp/markdup -s -f ${meta.id}.markdup.stats \\
            --reference ${fasta} -O cram --write-index - ${meta.id}.md.cram
    # markdup writes <out>.crai with --write-index
    test -s ${meta.id}.md.cram.crai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa-mem2: \$(bwa-mem2 version 2>&1 | grep -v Looking | head -n1)
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
