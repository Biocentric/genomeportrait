//
// VCF_TRAITS: look up the curated non-clinical "popular traits" SNP panel in the sample
//
include { TRAIT_LOOKUP } from '../../modules/local/trait_lookup'

workflow VCF_TRAITS {

    take:
    ch_vcf       // channel: [ meta, vcf, tbi ]   (annotated VCF, carries rsIDs)
    trait_panel  // file: assets/trait_snps.tsv

    main:
    ch_versions = Channel.empty()

    TRAIT_LOOKUP ( ch_vcf, trait_panel )
    ch_versions = ch_versions.mix(TRAIT_LOOKUP.out.versions.first())

    emit:
    results  = TRAIT_LOOKUP.out.tsv
    versions = ch_versions
}
