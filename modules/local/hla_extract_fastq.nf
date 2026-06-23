process HLA_EXTRACT_FASTQ {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::samtools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(cram), path(crai)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("${meta.id}_hla_1.fq.gz"), path("${meta.id}_hla_2.fq.gz"), emit: reads
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // T1K bundles samtools 0.1.19 (no collate/fastq), so extract the MHC reads to paired
    // FASTQ here with a modern samtools and hand T1K FASTQ input (needs only -f).
    // ALT-aware aligners (bwa-mem2) route divergent HLA reads to the dedicated HLA-* and
    // chr6_*_alt decoy contigs, so collect those in addition to the primary MHC interval.
    // The MHC contig name varies by reference convention (chr6 vs 6); detect from the header.
    """
    samtools view -H -T $fasta $cram > hdr.sam
    MHC=\$(awk -F'\\t' '/^@SQ/{for(i=1;i<=NF;i++){if(\$i=="SN:chr6"){print "chr6"} else if(\$i=="SN:6"){print "6"}}}' hdr.sam | head -1)
    [ -n "\$MHC" ] || MHC=chr6
    # primary MHC interval + every HLA-* and chr6_*_alt decoy contig present in the header
    ALTS=\$(sed -n 's/.*\\tSN:\\(HLA-[^\\t]*\\).*/\\1/p; s/.*\\tSN:\\(chr6_[A-Za-z0-9]*_alt\\).*/\\1/p' hdr.sam | tr '\\n' ' ')
    samtools view -@ $task.cpus -b -T $fasta $cram \${MHC}:28510120-33480577 \$ALTS > mhc.bam
    samtools collate -@ $task.cpus -u -O mhc.bam | \\
        samtools fastq -@ $task.cpus \\
            -1 ${meta.id}_hla_1.fq.gz -2 ${meta.id}_hla_2.fq.gz \\
            -0 /dev/null -s /dev/null -n

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
