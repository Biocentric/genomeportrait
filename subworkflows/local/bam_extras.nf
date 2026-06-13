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

    // Combine per sample, keyed on the string id (not the meta map) so absent
    // telomere/mito results become empty-file placeholders rather than corrupting
    // the tuple. Missing values -> [] (no staged file); the report script guards on them.
    ch_report_in = ch_bam.map { meta, b, i -> [ meta.id, meta ] }
        .join( ch_telseq.map { meta, f -> [ meta.id, f ] }, remainder: true )
        .join( ch_mito.map   { meta, f -> [ meta.id, f ] }, remainder: true )
        .map  { id, meta, tel, mito -> [ meta, tel ?: [], mito ?: [] ] }

    EXTRAS_REPORT ( ch_report_in )
    ch_versions = ch_versions.mix(EXTRAS_REPORT.out.versions.first())

    emit:
    results  = EXTRAS_REPORT.out.tsv
    versions = ch_versions
}
