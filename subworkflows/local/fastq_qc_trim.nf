//
// FASTQ_QC_TRIM: adapter/quality trimming (fastp) + FastQC
//
include { FASTP  } from '../../modules/local/fastp'
include { FASTQC } from '../../modules/local/fastqc'

workflow FASTQ_QC_TRIM {

    take:
    ch_reads      // channel: [ meta, [ reads ] ]
    skip_trimming // value

    main:
    ch_versions = Channel.empty()
    ch_multiqc  = Channel.empty()

    FASTQC ( ch_reads )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    ch_multiqc  = ch_multiqc.mix(FASTQC.out.zip)

    if (!skip_trimming) {
        FASTP ( ch_reads )
        ch_trimmed  = FASTP.out.reads
        ch_versions = ch_versions.mix(FASTP.out.versions.first())
        ch_multiqc  = ch_multiqc.mix(FASTP.out.json)
    } else {
        ch_trimmed = ch_reads
    }

    emit:
    reads    = ch_trimmed
    multiqc  = ch_multiqc
    versions = ch_versions
}
