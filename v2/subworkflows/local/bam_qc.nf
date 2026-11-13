//
// BAM_QC: samtools stats + Qualimap. Coverage (mosdepth) runs in the main workflow, because
// GENOME_FACTS needs it whether or not QC is skipped.
//
include { SAMTOOLS_STATS } from '../../modules/local/samtools_stats'
include { QUALIMAP_BAMQC } from '../../modules/local/qualimap_bamqc'

workflow BAM_QC {

    take:
    ch_bam          // channel: [ meta, bam/cram, index ]  analysis alignment
    ch_qualimap_bam // channel: [ meta, bam, bai ]         a real BAM; samples without one skip Qualimap
    ch_reference    // map

    main:
    ch_versions = Channel.empty()

    SAMTOOLS_STATS ( ch_bam, ch_reference.fasta )

    // Qualimap 2.3 cannot read CRAM ("BAM file is empty or corrupt"), so it runs on a BAM
    // that already exists (MarkDuplicates or Parabricks output, or a samplesheet BAM) and is
    // skipped for samples that only have a CRAM. Converting 40x back to BAM for one QC
    // module is not worth ~80 GB; mosdepth and samtools stats cover the same ground.
    QUALIMAP_BAMQC ( ch_qualimap_bam )

    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first(), QUALIMAP_BAMQC.out.versions.first())

    emit:
    multiqc  = SAMTOOLS_STATS.out.stats.mix(QUALIMAP_BAMQC.out.results)
    versions = ch_versions
}
