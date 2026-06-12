//
// BAM_FORENSIC_STR: CODIS core STR profile (HipSTR) + pathogenic-locus repeat sizing (ExpansionHunter)
//   The CODIS profile is the same set of markers used in forensic identity databases.
//
include { HIPSTR          } from '../../modules/local/hipstr'
include { EXPANSIONHUNTER } from '../../modules/local/expansionhunter'
include { STR_REPORT      } from '../../modules/local/str_report'

workflow BAM_FORENSIC_STR {

    take:
    ch_bam       // channel: [ meta, cram, crai ]
    ch_reference // map

    main:
    ch_versions = Channel.empty()

    HIPSTR ( ch_bam, ch_reference.fasta, ch_reference.fai, ch_reference.hipstr_codis )
    EXPANSIONHUNTER ( ch_bam, ch_reference.fasta, ch_reference.fai, ch_reference.str_catalog )
    ch_versions = ch_versions.mix(HIPSTR.out.versions.first(), EXPANSIONHUNTER.out.versions.first())

    STR_REPORT ( HIPSTR.out.vcf.join(EXPANSIONHUNTER.out.vcf) )
    ch_versions = ch_versions.mix(STR_REPORT.out.versions.first())

    emit:
    results  = STR_REPORT.out.tsv
    versions = ch_versions
}
