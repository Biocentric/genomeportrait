//
// VCF_ANCESTRY: PCA projection onto 1000G+HGDP, supervised ADMIXTURE, ROH/inbreeding, KING kinship
//
include { PLINK2_MERGE_REF      } from '../../modules/local/plink2_merge_ref'
include { PLINK2_PCA_PROJECT    } from '../../modules/local/plink2_pca_project'
include { ADMIXTURE_SUPERVISED  } from '../../modules/local/admixture_supervised'
include { PLINK2_ROH            } from '../../modules/local/plink2_roh'
include { KING_KINSHIP          } from '../../modules/local/king_kinship'
include { ANCESTRY_REPORT       } from '../../modules/local/ancestry_report'

workflow VCF_ANCESTRY {

    take:
    ch_vcf       // channel: [ meta, vcf, tbi ]
    ch_reference // map (must contain kgp_hgdp_panel)

    main:
    ch_versions = Channel.empty()

    // Merge the sample's genotypes with the reference panel on shared sites
    PLINK2_MERGE_REF ( ch_vcf, ch_reference.kgp_hgdp_panel )
    ch_versions = ch_versions.mix(PLINK2_MERGE_REF.out.versions.first())

    // Project the sample into reference PCA space
    PLINK2_PCA_PROJECT ( PLINK2_MERGE_REF.out.merged, ch_reference.kgp_hgdp_panel )
    ch_versions = ch_versions.mix(PLINK2_PCA_PROJECT.out.versions.first())

    // Supervised global ancestry proportions
    ADMIXTURE_SUPERVISED ( PLINK2_MERGE_REF.out.merged, ch_reference.kgp_hgdp_panel )
    ch_versions = ch_versions.mix(ADMIXTURE_SUPERVISED.out.versions.first())

    // Runs of homozygosity / inbreeding coefficient F
    PLINK2_ROH ( ch_vcf )
    ch_versions = ch_versions.mix(PLINK2_ROH.out.versions.first())

    // KING kinship vs reference (and self-QC)
    KING_KINSHIP ( PLINK2_MERGE_REF.out.merged )
    ch_versions = ch_versions.mix(KING_KINSHIP.out.versions.first())

    // Assemble ancestry tables + PCA/admixture plots
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
