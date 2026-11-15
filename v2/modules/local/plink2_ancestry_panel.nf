process PLINK2_ANCESTRY_PANEL {
    tag "ancestry_snps"
    label 'process_medium'
    label 'process_long'
    // Param-keyed storeDir: rebuilds if the SNP-selection knobs change, else reused.
    // "hgdp1kgp_v2": the panel and the sites format changed; an older cache must not be reused
    storeDir "${params.reference_base}/ancestry_snps/hgdp1kgp_v2_maf${params.ancestry_min_maf}_n${params.ancestry_max_snps}"

    conda "bioconda::plink2=2.00a5.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0' :
        'biocontainers/plink2:2.00a5.10--h4ac6f70_0' }"

    input:
    tuple path(ref_pgen), path(ref_pvar), path(ref_psam), path(related)

    output:
    tuple path("panel_anc.bed"), path("panel_anc.bim"), path("panel_anc.fam"), emit: panel
    path "anc_sites.tsv",                                                       emit: sites
    path "versions.yml",                                                        emit: versions

    script:
    // .pvar.zst (the HGDP+1kGP release) needs the vzs modifier; relatives are dropped, as for PGS.
    def vzs = ref_pvar.name.endsWith('.zst') ? 'vzs' : ''
    // Define the ancestry SNP set ONCE from the panel (sample-independent, cached):
    // autosomal biallelic common SNPs, allele-aware IDs (@:#:ref:alt), LD-pruned, capped.
    // Emit BED so plink 1.9 can later --bmerge it. anc_sites.tsv = chr-PREFIXED CHROM POS REF,ALT
    // (+ alleles) for genotyping the sample BAM.
    """
    plink2 --pfile ${ref_pgen.baseName} ${vzs} --remove ${related} --allow-extra-chr --autosome \\
        --snps-only --max-alleles 2 --min-alleles 2 \\
        --maf ${params.ancestry_min_maf} \\
        --set-all-var-ids '@:#:\$r:\$a' --new-id-max-allele-len 100 --rm-dup force-first \\
        --indep-pairwise 200 50 0.2 --out prune

    # Every k-th pruned SNP, so the set spans all autosomes. (v1 took the first N lines of
    # prune.in, which is in genome order: chromosomes 1-3 only.)
    total=\$(wc -l < prune.prune.in)
    step=\$(( (total + ${params.ancestry_max_snps} - 1) / ${params.ancestry_max_snps} ))
    awk -v k=\$step 'NR % k == 1 || k == 1' prune.prune.in > keep.ids
    echo "ancestry SNPs: \$(wc -l < keep.ids) of \$total pruned (every \${step})" >&2
    # SNP-only here too: otherwise long indel allele codes blow past --new-id-max-allele-len
    # before --extract can drop them.
    plink2 --pfile ${ref_pgen.baseName} ${vzs} --remove ${related} --allow-extra-chr \\
        --snps-only --max-alleles 2 --min-alleles 2 \\
        --set-all-var-ids '@:#:\$r:\$a' --new-id-max-allele-len 100 --rm-dup force-first \\
        --extract keep.ids --make-bed --out panel_anc

    # sites for bcftools mpileup (.bim: CHROM ID cM POS A1 A2 -> chr-prefixed CHROM POS A1 A2)
    # CHROM POS REF,ALT: the alleles come from the variant ID (chr:pos:ref:alt), so the sample
    # can be genotyped for exactly these alleles (bcftools call -C alleles) and hom-ref sites
    # get IDs that match the panel.
    awk '{ split(\$2, a, ":"); print "chr"\$1"\\t"\$4"\\t"a[3]","a[4] }' panel_anc.bim | sort -k1,1 -k2,2n > anc_sites.tsv
    echo "ancestry SNPs: \$(wc -l < anc_sites.tsv)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    touch "panel_anc.bed"
    touch "panel_anc.bim"
    touch "panel_anc.fam"
    touch "anc_sites.tsv"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: true
    END_VERSIONS
    """
}
