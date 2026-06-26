//
// BAM_ANCESTRY: genotype the sample at a curated common-SNP ancestry panel directly from the
// alignment, merge into the 1000G+HGDP reference, then PCA projection, supervised ADMIXTURE,
// ROH/inbreeding and KING kinship. (Genotyping from the BAM — not the variant-only VCF — and a
// MAF-filtered SNP set are both required for a single 8x sample to land in the right cluster.)
//
include { PLINK2_ANCESTRY_PANEL   } from '../../modules/local/plink2_ancestry_panel'
include { BCFTOOLS_GENOTYPE_SITES } from '../../modules/local/bcftools_genotype_sites'
include { PLINK2_ANCESTRY_SAMPLE  } from '../../modules/local/plink2_ancestry_sample'
include { PLINK_ANCESTRY_MERGE    } from '../../modules/local/plink_ancestry_merge'
include { PLINK2_PCA_PROJECT      } from '../../modules/local/plink2_pca_project'
include { ADMIXTURE_SUPERVISED    } from '../../modules/local/admixture_supervised'
include { PLINK2_ROH              } from '../../modules/local/plink2_roh'
include { KING_KINSHIP            } from '../../modules/local/king_kinship'
include { ANCESTRY_REPORT         } from '../../modules/local/ancestry_report'

workflow BAM_ANCESTRY {

    take:
    ch_bam       // channel: [ meta, cram, crai ]
    ch_vcf       // channel: [ meta, vcf, tbi ]   (for ROH)
    ch_reference // map (must contain kgp_hgdp_panel = [pgen,pvar,psam])

    main:
    ch_versions = Channel.empty()

    // 1. Define the ancestry SNP set from the panel (common, LD-pruned, allele-aware; cached)
    PLINK2_ANCESTRY_PANEL ( ch_reference.kgp_hgdp_panel )
    ch_panel = PLINK2_ANCESTRY_PANEL.out.panel
    ch_versions = ch_versions.mix(PLINK2_ANCESTRY_PANEL.out.versions)

    // 2. Genotype the sample at those sites from the alignment (hom-ref captured)
    BCFTOOLS_GENOTYPE_SITES ( ch_bam, ch_reference.fasta, ch_reference.fai, PLINK2_ANCESTRY_PANEL.out.sites )
    ch_versions = ch_versions.mix(BCFTOOLS_GENOTYPE_SITES.out.versions.first())

    // 3. Sample -> bed at the panel SNPs, then 4. merge into the reference panel
    PLINK2_ANCESTRY_SAMPLE ( BCFTOOLS_GENOTYPE_SITES.out.vcf, ch_panel )
    PLINK_ANCESTRY_MERGE   ( PLINK2_ANCESTRY_SAMPLE.out.bed, ch_panel )
    ch_comb = PLINK_ANCESTRY_MERGE.out.combined
    ch_versions = ch_versions.mix(PLINK2_ANCESTRY_SAMPLE.out.versions.first(), PLINK_ANCESTRY_MERGE.out.versions.first())

    // 5. PCA projection onto the reference basis
    PLINK2_PCA_PROJECT ( ch_comb, ch_panel )
    ch_versions = ch_versions.mix(PLINK2_PCA_PROJECT.out.versions.first())

    // 6. Supervised global ancestry proportions
    ADMIXTURE_SUPERVISED ( ch_comb, ch_panel, ch_reference.kgp_hgdp_panel.map { it[2] } )
    ch_versions = ch_versions.mix(ADMIXTURE_SUPERVISED.out.versions.first())

    // 7. Runs of homozygosity / inbreeding F (from the sample's own variants)
    PLINK2_ROH ( ch_vcf )
    ch_versions = ch_versions.mix(PLINK2_ROH.out.versions.first())

    // 8. KING kinship vs the reference
    KING_KINSHIP ( ch_comb )
    ch_versions = ch_versions.mix(KING_KINSHIP.out.versions.first())

    // 9. Assemble ancestry tables + PCA/admixture plots
    ANCESTRY_REPORT (
        PLINK2_PCA_PROJECT.out.eigenvec
            .join(ADMIXTURE_SUPERVISED.out.q)
            .join(PLINK2_ROH.out.roh)
            .join(KING_KINSHIP.out.kin),
        ch_reference.kgp_hgdp_panel
    )
    ch_versions = ch_versions.mix(ANCESTRY_REPORT.out.versions.first())

    emit:
    results  = ANCESTRY_REPORT.out.bundle
    versions = ch_versions
}
