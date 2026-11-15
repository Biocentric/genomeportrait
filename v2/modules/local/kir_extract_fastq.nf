process KIR_EXTRACT_FASTQ {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::samtools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(cram), path(crai)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("${meta.id}_kir_1.fq.gz"), path("${meta.id}_kir_2.fq.gz"), emit: reads
    path "versions.yml",                                                           emit: versions

    script:
    // The KIR cluster (leukocyte receptor complex, chr19:54.7-54.9 Mb in GRCh38) plus every
    // chr19 ALT contig: GRCh38 carries many alternative KIR haplotypes, and an alignment that
    // is not ALT-aware puts KIR reads on them. T1K re-aligns everything to its own allele
    // database, so it does not matter where the reads were placed, only that they are passed.
    """
    samtools view -H -T ${fasta} ${cram} > hdr.sam
    LRC=\$(awk -F'\\t' '/^@SQ/{ for (i = 1; i <= NF; i++) { if (\$i == "SN:chr19") print "chr19"; else if (\$i == "SN:19") print "19" } }' hdr.sam | head -1)
    [ -n "\$LRC" ] || LRC=chr19
    ALTS=\$(sed -n 's/.*\\tSN:\\(chr19_[A-Za-z0-9]*_alt\\).*/\\1/p' hdr.sam | tr '\\n' ' ')
    samtools view -@ ${task.cpus} -b -T ${fasta} ${cram} \${LRC}:54700000-54900000 \$ALTS > kir.bam
    samtools collate -@ ${task.cpus} -u -O kir.bam collate_tmp | \\
        samtools fastq -@ ${task.cpus} -1 ${meta.id}_kir_1.fq.gz -2 ${meta.id}_kir_2.fq.gz -0 /dev/null -s /dev/null -n
    rm -f kir.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
