process PLINK2_ANCESTRY_SAMPLE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::plink2=2.00a5.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0' :
        'biocontainers/plink2:2.00a5.10--h4ac6f70_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple path(panel_bed), path(panel_bim), path(panel_fam)

    output:
    tuple val(meta), path("${meta.id}.sample.bed"), path("${meta.id}.sample.bim"), path("${meta.id}.sample.fam"), emit: bed
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Sample (BAM-genotyped) VCF -> bed with the SAME allele-aware IDs as the panel
    // (plink2 normalises chr1->1 so they line up), keeping only the panel's ancestry SNPs.
    """
    awk '{print \$2}' $panel_bim > panel.ids
    plink2 --vcf $vcf --allow-extra-chr --snps-only --max-alleles 2 --min-alleles 2 \\
        --set-all-var-ids '@:#:\$r:\$a' --new-id-max-allele-len 100 --rm-dup force-first \\
        --extract panel.ids --make-bed --out ${meta.id}.sample

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ .*//')
    END_VERSIONS
    """
}
