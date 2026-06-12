//
// BAM_EXTRAS: telomere length estimate (Telseq) + mtDNA heteroplasmy (mutserve)
//
include { TELSEQ   } from '../../modules/local/telseq'
include { MUTSERVE } from '../../modules/local/mutserve'
include { EXTRAS_REPORT } from '../../modules/local/extras_report'

workflow BAM_EXTRAS {

    take:
    ch_bam       // channel: [ meta, cram, crai ]
    ch_reference // map

    main:
    ch_versions = Channel.empty()

    ch_telseq = Channel.empty()
    ch_mito   = Channel.empty()

    if (!params.skip_telomere) {
        TELSEQ ( ch_bam )
        ch_telseq   = TELSEQ.out.tsv
        ch_versions = ch_versions.mix(TELSEQ.out.versions.first())
    }
    if (!params.skip_mito_heteroplasmy) {
        MUTSERVE ( ch_bam, ch_reference.fasta, ch_reference.fai )
        ch_mito     = MUTSERVE.out.txt
        ch_versions = ch_versions.mix(MUTSERVE.out.versions.first())
    }

    EXTRAS_REPORT (
        ch_bam.map{ meta, b, i -> [meta] }
            .join(ch_telseq, remainder: true)
            .join(ch_mito, remainder: true)
    )
    ch_versions = ch_versions.mix(EXTRAS_REPORT.out.versions.first())

    emit:
    results  = EXTRAS_REPORT.out.tsv
    versions = ch_versions
}
