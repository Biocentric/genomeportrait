//
// BAM_QC: coverage (mosdepth), samtools stats, and optionally Qualimap
//
include { MOSDEPTH                      } from '../../modules/local/mosdepth'
include { SAMTOOLS_STATS                } from '../../modules/local/samtools_stats'
include { CRAM_TO_BAM as QC_CRAM_TO_BAM } from '../../modules/local/cram_to_bam'
include { QUALIMAP_BAMQC                } from '../../modules/local/qualimap_bamqc'

workflow BAM_QC {

    take:
    ch_bam       // channel: [ meta, bam/cram, bai/crai ]
    ch_reference // map

    main:
    ch_versions = Channel.empty()
    ch_multiqc  = Channel.empty()

    MOSDEPTH ( ch_bam, ch_reference.fasta )
    SAMTOOLS_STATS ( ch_bam, ch_reference.fasta )
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions.first(), SAMTOOLS_STATS.out.versions.first())
    ch_multiqc  = ch_multiqc
        .mix(MOSDEPTH.out.summary)
        .mix(MOSDEPTH.out.global_dist)
        .mix(SAMTOOLS_STATS.out.stats)

    // Qualimap 2.3 cannot read CRAM (it reports "BAM file is empty or corrupt"), so it needs
    // a decompressed BAM — roughly 3-4x the CRAM size, i.e. >100 GB at production depth.
    // mosdepth + samtools stats already cover coverage/insert-size/GC/MAPQ in MultiQC, so
    // this is opt-in: enable with --skip_qualimap false when the disk can take it.
    if (!params.skip_qualimap) {
        QC_CRAM_TO_BAM ( ch_bam, ch_reference.fasta, ch_reference.fai )
        QUALIMAP_BAMQC ( QC_CRAM_TO_BAM.out.bam )
        ch_versions = ch_versions.mix(QC_CRAM_TO_BAM.out.versions.first(), QUALIMAP_BAMQC.out.versions.first())
        ch_multiqc  = ch_multiqc.mix(QUALIMAP_BAMQC.out.results)
    }

    emit:
    multiqc  = ch_multiqc
    versions = ch_versions
}
