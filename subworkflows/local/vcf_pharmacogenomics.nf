//
// BAM_PHARMACOGENOMICS: PharmCAT star-allele calling + drug-gene annotations
//   Reports metaboliser phenotypes (e.g. CYP2D6, CYP2C19) and CPIC guidance.
//   Informational only — not a prescription tool.
//
include { PHARMCAT        } from '../../modules/local/pharmcat'
include { PHARMCAT_REPORT } from '../../modules/local/pharmcat_report'

workflow BAM_PHARMACOGENOMICS {

    take:
    ch_bam       // channel: [ meta, bam/cram, bai/crai ]
    ch_reference // map

    main:
    ch_versions = Channel.empty()

    PHARMCAT ( ch_bam, ch_reference.fasta, ch_reference.fai )
    ch_versions = ch_versions.mix(PHARMCAT.out.versions.first())

    PHARMCAT_REPORT ( PHARMCAT.out.json )
    ch_versions = ch_versions.mix(PHARMCAT_REPORT.out.versions.first())

    emit:
    results  = PHARMCAT_REPORT.out.tsv
    versions = ch_versions
}
