//
// BAM_QC: coverage (mosdepth), samtools stats, Qualimap
//
include { MOSDEPTH       } from '../../modules/local/mosdepth'
include { SAMTOOLS_STATS } from '../../modules/local/samtools_stats'
include { QUALIMAP_BAMQC } from '../../modules/local/qualimap_bamqc'

workflow BAM_QC {

    take:
    ch_bam       // channel: [ meta, bam/cram, bai/crai ]  (analysis alignment, CRAM after BQSR)
    ch_qualimap_bam // channel: [ meta, bam, bai ]  — a real BAM (MarkDuplicates output)
    ch_reference // map

    main:
    ch_versions = Channel.empty()
    ch_multiqc  = Channel.empty()

    MOSDEPTH ( ch_bam, ch_reference.fasta )
    SAMTOOLS_STATS ( ch_bam, ch_reference.fasta )

    // Qualimap 2.3 cannot read CRAM ("BAM file is empty or corrupt"), so run it on the
    // MarkDuplicates BAM, which already exists in the work dir — no extra conversion or
    // disk. BQSR only rewrites base-quality scores, so coverage/insert-size/GC/MAPQ
    // metrics are the same as on the recalibrated alignment.
    QUALIMAP_BAMQC ( ch_qualimap_bam )

    ch_versions = ch_versions.mix(MOSDEPTH.out.versions.first(), SAMTOOLS_STATS.out.versions.first(), QUALIMAP_BAMQC.out.versions.first())
    ch_multiqc  = ch_multiqc
        .mix(MOSDEPTH.out.summary)
        .mix(MOSDEPTH.out.global_dist)
        .mix(SAMTOOLS_STATS.out.stats)
        .mix(QUALIMAP_BAMQC.out.results)

    emit:
    multiqc  = ch_multiqc
    versions = ch_versions
}
