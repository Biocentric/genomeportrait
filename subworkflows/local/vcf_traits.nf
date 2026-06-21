//
// VCF_TRAITS: genotype the curated non-clinical "popular traits" SNP panel directly from
// the alignment (targeted mpileup), so homozygous-reference sites are captured, then
// interpret each genotype against the panel.
//
include { TRAIT_GENOTYPE } from '../../modules/local/trait_genotype'
include { TRAIT_LOOKUP   } from '../../modules/local/trait_lookup'

workflow VCF_TRAITS {

    take:
    ch_bam       // channel: [ meta, bam/cram, bai/crai ]
    ch_reference // map
    trait_panel  // file: assets/trait_snps.tsv

    main:
    ch_versions = Channel.empty()

    TRAIT_GENOTYPE ( ch_bam, ch_reference.fasta, ch_reference.fai, trait_panel )
    TRAIT_LOOKUP   ( TRAIT_GENOTYPE.out.geno, trait_panel )
    ch_versions = ch_versions.mix(TRAIT_GENOTYPE.out.versions.first(), TRAIT_LOOKUP.out.versions.first())

    emit:
    results  = TRAIT_LOOKUP.out.tsv
    versions = ch_versions
}
