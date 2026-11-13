//
// BAM_MARKDUP_BQSR: GATK4 MarkDuplicates, then either BQSR (opt-in) or a plain CRAM conversion.
// Alignments already duplicate-marked in the streamed aligner (BWAMEM2_MARKDUP_CRAM) join
// after the MarkDuplicates step: they get BQSR when it is on and are passed through otherwise.
//
// BQSR is off by default in v2: DeepVariant is trained on unrecalibrated base qualities and
// gains nothing from it, and it costs one more whole-genome pass. Either way the analysis
// alignment that leaves this subworkflow is a CRAM.
//
include { GATK4_MARKDUPLICATES } from '../../modules/local/gatk4_markduplicates'
include { GATK4_BASERECALIBRATOR } from '../../modules/local/gatk4_baserecalibrator'
include { GATK4_APPLYBQSR       } from '../../modules/local/gatk4_applybqsr'
include { SAMTOOLS_BAM2CRAM     } from '../../modules/local/samtools_bam2cram'

workflow BAM_MARKDUP_BQSR {

    take:
    ch_bam       // channel: [ meta, bam, bai ]     sorted, duplicates not yet marked
    ch_reference // map
    skip_bqsr    // value
    ch_marked    // channel: [ meta, cram, crai ]   already duplicate-marked (streamed path)

    main:
    ch_versions = Channel.empty()

    GATK4_MARKDUPLICATES ( ch_bam, ch_reference.fasta, ch_reference.fai )
    ch_md       = GATK4_MARKDUPLICATES.out.bam
    ch_metrics  = GATK4_MARKDUPLICATES.out.metrics
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions.first())

    // GATK reads CRAM with the reference, so streamed CRAMs go through BQSR unchanged.
    ch_to_bqsr = ch_md.mix(ch_marked)

    if (!skip_bqsr) {
        GATK4_BASERECALIBRATOR (
            ch_to_bqsr,
            ch_reference.fasta, ch_reference.fai, ch_reference.dict,
            ch_reference.dbsnp, ch_reference.known_indels, ch_reference.mills
        )
        GATK4_APPLYBQSR (
            ch_to_bqsr.join(GATK4_BASERECALIBRATOR.out.table),
            ch_reference.fasta, ch_reference.fai, ch_reference.dict
        )
        ch_final    = GATK4_APPLYBQSR.out.cram
        ch_versions = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions.first(), GATK4_APPLYBQSR.out.versions.first())
    } else {
        SAMTOOLS_BAM2CRAM ( ch_md, ch_reference.fasta, ch_reference.fai )
        ch_final    = SAMTOOLS_BAM2CRAM.out.cram.mix(ch_marked)
        ch_versions = ch_versions.mix(SAMTOOLS_BAM2CRAM.out.versions.first())
    }

    emit:
    bam      = ch_final   // [ meta, cram, crai ]
    md_bam   = ch_md      // [ meta, bam, bai ] — real BAM, for tools that can't read CRAM
    metrics  = ch_metrics
    versions = ch_versions
}
