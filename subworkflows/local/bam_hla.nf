//
// BAM_HLA: HLA class I & II typing at 4-digit (G-group) resolution with T1K
//
include { HLA_T1K     } from '../../modules/local/hla_t1k'
include { HLA_REPORT  } from '../../modules/local/hla_report'

workflow BAM_HLA {

    take:
    ch_bam       // channel: [ meta, cram, crai ]
    ch_reference // map

    main:
    ch_versions = Channel.empty()

    HLA_T1K ( ch_bam, ch_reference.fasta, ch_reference.fai )
    ch_versions = ch_versions.mix(HLA_T1K.out.versions.first())

    HLA_REPORT ( HLA_T1K.out.alleles )
    ch_versions = ch_versions.mix(HLA_REPORT.out.versions.first())

    emit:
    results  = HLA_REPORT.out.tsv
    versions = ch_versions
}
