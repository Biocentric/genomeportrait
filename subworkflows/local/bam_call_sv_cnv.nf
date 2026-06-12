//
// BAM_CALL_SV_CNV: structural variants (Manta) + copy-number (CNVkit)
//
include { MANTA_GERMLINE } from '../../modules/local/manta_germline'
include { CNVKIT_BATCH   } from '../../modules/local/cnvkit_batch'
include { SVCNV_SUMMARY  } from '../../modules/local/svcnv_summary'

workflow BAM_CALL_SV_CNV {

    take:
    ch_bam       // channel: [ meta, cram, crai ]
    ch_reference // map

    main:
    ch_versions = Channel.empty()
    ch_sv  = Channel.empty()
    ch_cnv = Channel.empty()

    if (params.call_sv) {
        MANTA_GERMLINE ( ch_bam, ch_reference.fasta, ch_reference.fai )
        ch_sv       = MANTA_GERMLINE.out.sv_vcf
        ch_versions = ch_versions.mix(MANTA_GERMLINE.out.versions.first())
    }
    if (params.call_cnv) {
        CNVKIT_BATCH ( ch_bam, ch_reference.fasta )
        ch_cnv      = CNVKIT_BATCH.out.cnr
        ch_versions = ch_versions.mix(CNVKIT_BATCH.out.versions.first())
    }

    // Join SV + CNV per sample (either may be empty) and build a compact summary table
    SVCNV_SUMMARY (
        ch_bam.map{ meta, bam, bai -> [meta] }
            .join(ch_sv.map{ m,f,t -> [m,f] }, remainder: true)
            .join(ch_cnv, remainder: true)
    )
    ch_versions = ch_versions.mix(SVCNV_SUMMARY.out.versions.first())

    emit:
    sv       = ch_sv
    cnv      = ch_cnv
    summary  = SVCNV_SUMMARY.out.tsv
    versions = ch_versions
}
