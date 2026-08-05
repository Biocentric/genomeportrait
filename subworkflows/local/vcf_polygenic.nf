//
// VCF_POLYGENIC: polygenic scores from the PGS Catalog
//   Engine mirrors pgsc_calc: harmonise scoring files to the target build, then
//   apply additive weights with PLINK2 --score. Percentiles are derived against
//   the 1000G distribution shipped with each scoring file where available.
//
include { PGSCATALOG_DOWNLOAD } from '../../modules/local/pgscatalog_download'
include { PLINK2_VCF2PGEN     } from '../../modules/local/plink2_vcf2pgen'
include { PLINK2_SCORE        } from '../../modules/local/plink2_score'
include { PLINK2_SCORE_REFERENCE } from '../../modules/local/plink2_score_reference'
include { PRS_REPORT          } from '../../modules/local/prs_report'

workflow VCF_POLYGENIC {

    take:
    ch_vcf       // channel: [ meta, vcf, tbi ]
    ch_reference // map

    main:
    ch_versions = Channel.empty()

    // Fetch + harmonise the requested PGS Catalog scoring files (one combined file)
    PGSCATALOG_DOWNLOAD ( params.pgs_ids )
    ch_versions = ch_versions.mix(PGSCATALOG_DOWNLOAD.out.versions)

    PLINK2_VCF2PGEN ( ch_vcf )
    PLINK2_SCORE ( PLINK2_VCF2PGEN.out.pgen, PGSCATALOG_DOWNLOAD.out.scorefiles )
    ch_versions = ch_versions.mix(PLINK2_VCF2PGEN.out.versions.first(), PLINK2_SCORE.out.versions.first())

    // Reference distribution (cached, sample-independent) so raw scores become percentiles
    ch_ref_scores = Channel.value([])
    if (!params.skip_ancestry) {
        PLINK2_SCORE_REFERENCE ( ch_reference.kgp_hgdp_panel, PGSCATALOG_DOWNLOAD.out.scorefiles )
        ch_ref_scores = PLINK2_SCORE_REFERENCE.out.scores
        ch_versions   = ch_versions.mix(PLINK2_SCORE_REFERENCE.out.versions)
    }

    PRS_REPORT ( PLINK2_SCORE.out.sscore, PGSCATALOG_DOWNLOAD.out.metadata, ch_ref_scores )
    ch_versions = ch_versions.mix(PRS_REPORT.out.versions.first())

    emit:
    scores   = PRS_REPORT.out.tsv
    versions = ch_versions
}
