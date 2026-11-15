process RHD_DEPTH {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::samtools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fai
    path  regions     // assets/blood_groups/rhd_regions.bed

    output:
    tuple val(meta), path("${meta.id}.rhd_depth.tsv"), emit: tsv
    path "versions.yml",                               emit: versions

    script:
    // RhD-negative is mostly a deletion of the whole RHD gene. Reads per base over RHD against
    // two flanking control windows: ~1 with two copies, ~0.5 with one, ~0 with none. MAPQ >= 20
    // keeps reads from the paralogous RHCE out; duplicates, secondary and QC-failed reads are
    // not counted.
    """
    : > counts.tsv
    while IFS=\$'\\t' read -r chrom start end name; do
        n=\$(samtools view -c -q 20 -F 0xD04 -T ${fasta} ${bam} "\${chrom}:\$((start + 1))-\${end}")
        printf "%s\\t%s\\t%s\\n" "\$name" "\$n" "\$((end - start))" >> counts.tsv
    done < ${regions}
    awk -F'\\t' '\$1 == "RHD" { r = \$2 / \$3 } \$1 ~ /^ctrl/ { c += \$2; l += \$3 }
        END { printf "metric\\tvalue\\nrhd_reads_per_kb\\t%.3f\\ncontrol_reads_per_kb\\t%.3f\\nratio\\t%.4f\\n", 1000 * r, 1000 * c / l, (c > 0 ? r / (c / l) : 0) }' \\
        counts.tsv > ${meta.id}.rhd_depth.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
