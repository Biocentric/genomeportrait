//
// BAM_CALL_SV_CNV: structural variants (Manta) + copy-number (CNVpytor)
//   CNVpytor (read-depth, CNVnator successor) is efficient for single-sample germline
//   WGS, unlike CNVkit which targets exome / panel-of-normals designs.
//
include { MANTA_GERMLINE } from '../../modules/local/manta_germline'
include { CNVPYTOR       } from '../../modules/local/cnvpytor'
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
        CNVPYTOR ( ch_bam, ch_reference.fasta, ch_reference.fai )
        ch_cnv      = CNVPYTOR.out.calls
        ch_versions = ch_versions.mix(CNVPYTOR.out.versions.first())
    }

    // Join SV + CNV per sample (either may be empty) keyed on the string id, so an
    // absent caller becomes an empty-file placeholder instead of corrupting the tuple.
    ch_summary_in = ch_bam.map { meta, bam, bai -> [ meta.id, meta ] }
        .join( ch_sv.map  { m, f, t -> [ m.id, f ] }, remainder: true )
        .join( ch_cnv.map { m, f    -> [ m.id, f ] }, remainder: true )
        .map  { id, meta, sv, cnv -> [ meta, sv ?: [], cnv ?: [] ] }

    SVCNV_SUMMARY ( ch_summary_in )
    ch_versions = ch_versions.mix(SVCNV_SUMMARY.out.versions.first())

    emit:
    sv       = ch_sv
    cnv      = ch_cnv
    summary  = SVCNV_SUMMARY.out.tsv
    versions = ch_versions
}
