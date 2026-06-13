//
// BAM_MARKDUP_BQSR: GATK4 MarkDuplicates + Base Quality Score Recalibration
//
include { GATK4_MARKDUPLICATES } from '../../modules/local/gatk4_markduplicates'
include { GATK4_BASERECALIBRATOR } from '../../modules/local/gatk4_baserecalibrator'
include { GATK4_APPLYBQSR       } from '../../modules/local/gatk4_applybqsr'

workflow BAM_MARKDUP_BQSR {

    take:
    ch_bam       // channel: [ meta, bam, bai ]
    ch_reference // map
    skip_bqsr    // value

    main:
    ch_versions = Channel.empty()

    GATK4_MARKDUPLICATES ( ch_bam, ch_reference.fasta, ch_reference.fai )
    ch_md       = GATK4_MARKDUPLICATES.out.bam
    ch_metrics  = GATK4_MARKDUPLICATES.out.metrics
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions.first())

    if (!skip_bqsr) {
        GATK4_BASERECALIBRATOR (
            ch_md,
            ch_reference.fasta, ch_reference.fai, ch_reference.dict,
            ch_reference.dbsnp, ch_reference.known_indels, ch_reference.mills
        )
        GATK4_APPLYBQSR (
            ch_md.join(GATK4_BASERECALIBRATOR.out.table),
            ch_reference.fasta, ch_reference.fai, ch_reference.dict
        )
        ch_final    = GATK4_APPLYBQSR.out.cram
        ch_versions = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions.first(), GATK4_APPLYBQSR.out.versions.first())
    } else {
        ch_final = ch_md
    }

    emit:
    bam      = ch_final   // [ meta, cram/bam, crai/bai ]
    metrics  = ch_metrics
    versions = ch_versions
}
