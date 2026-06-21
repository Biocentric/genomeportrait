//
// BAM_FORENSIC_STR: CODIS core STR profile (HipSTR) + pathogenic-locus repeat sizing (ExpansionHunter)
//   The CODIS profile is the same set of markers used in forensic identity databases.
//
include { CRAM_TO_BAM     } from '../../modules/local/cram_to_bam'
include { HIPSTR          } from '../../modules/local/hipstr'
include { EXPANSIONHUNTER } from '../../modules/local/expansionhunter'
include { STR_REPORT      } from '../../modules/local/str_report'

workflow BAM_FORENSIC_STR {

    take:
    ch_bam       // channel: [ meta, cram, crai ]
    ch_reference // map

    main:
    ch_versions = Channel.empty()

    // HipSTR (and ExpansionHunter) are happiest with BAM; the analysis file is a CRAM,
    // and HipSTR's older htslib reads 0 reads from it. Convert once and feed both.
    CRAM_TO_BAM ( ch_bam, ch_reference.fasta, ch_reference.fai )
    ch_strbam = CRAM_TO_BAM.out.bam
    ch_versions = ch_versions.mix(CRAM_TO_BAM.out.versions.first())

    HIPSTR ( ch_strbam, ch_reference.fasta, ch_reference.fai, ch_reference.hipstr_codis, ch_reference.hipstr_ready )
    EXPANSIONHUNTER ( ch_strbam, ch_reference.fasta, ch_reference.fai, ch_reference.str_catalog )
    ch_versions = ch_versions.mix(HIPSTR.out.versions.first(), EXPANSIONHUNTER.out.versions.first())

    STR_REPORT ( HIPSTR.out.vcf.join(EXPANSIONHUNTER.out.vcf) )
    ch_versions = ch_versions.mix(STR_REPORT.out.versions.first())

    emit:
    results  = STR_REPORT.out.tsv
    versions = ch_versions
}
