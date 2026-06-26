process PLINK2_ANCESTRY_PANEL {
    tag "ancestry_snps"
    label 'process_medium'
    label 'process_long'
    // Param-keyed storeDir: rebuilds if the SNP-selection knobs change, else reused.
    storeDir "${params.reference_base}/ancestry_snps/maf${params.ancestry_min_maf}_n${params.ancestry_max_snps}"

    conda "bioconda::plink2=2.00a5.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0' :
        'biocontainers/plink2:2.00a5.10--h4ac6f70_0' }"

    input:
    tuple path(ref_pgen), path(ref_pvar), path(ref_psam)

    output:
    tuple path("panel_anc.bed"), path("panel_anc.bim"), path("panel_anc.fam"), emit: panel
    path "anc_sites.tsv",                                                       emit: sites
    path "versions.yml",                                                        emit: versions

    script:
    // Define the ancestry SNP set ONCE from the panel (sample-independent, cached):
    // autosomal biallelic common SNPs, allele-aware IDs (@:#:ref:alt), LD-pruned, capped.
    // Emit BED so plink 1.9 can later --bmerge it. anc_sites.tsv = chr-PREFIXED CHROM POS
    // (+ alleles) for genotyping the sample BAM.
    """
    plink2 --pfile ${ref_pgen.baseName} --allow-extra-chr --autosome \\
        --snps-only --max-alleles 2 --min-alleles 2 \\
        --maf ${params.ancestry_min_maf} \\
        --set-all-var-ids '@:#:\$r:\$a' --new-id-max-allele-len 100 --rm-dup force-first \\
        --indep-pairwise 200 50 0.2 --out prune

    head -n ${params.ancestry_max_snps} prune.prune.in > keep.ids
    plink2 --pfile ${ref_pgen.baseName} --allow-extra-chr \\
        --set-all-var-ids '@:#:\$r:\$a' --new-id-max-allele-len 100 \\
        --extract keep.ids --make-bed --out panel_anc

    # sites for bcftools mpileup (.bim: CHROM ID cM POS A1 A2 -> chr-prefixed CHROM POS A1 A2)
    awk '{print "chr"\$1"\\t"\$4"\\t"\$5"\\t"\$6}' panel_anc.bim | sort -k1,1 -k2,2n > anc_sites.tsv
    echo "ancestry SNPs: \$(wc -l < anc_sites.tsv)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ .*//')
    END_VERSIONS
    """
}
