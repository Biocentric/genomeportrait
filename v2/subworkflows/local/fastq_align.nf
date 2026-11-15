//
// FASTQ_ALIGN: align short reads (bwa-mem2) or long reads (minimap2); index.
//
// Both aligner modules pipe straight into `samtools sort`, so their output is already
// coordinate-sorted. v1 sorted it a second time here (SAMTOOLS_SORT): one redundant
// whole-genome pass and ~120 GB of scratch per 40x sample. v2 only indexes.
//
include { BWAMEM2_MEM    } from '../../modules/local/bwamem2_mem'
include { MINIMAP2_ALIGN } from '../../modules/local/minimap2_align'
include { SAMTOOLS_INDEX } from '../../modules/local/samtools_index'

workflow FASTQ_ALIGN {

    take:
    ch_reads     // channel: [ meta, [ reads ] ]
    ch_reference // map

    main:
    ch_versions = Channel.empty()

    // Branch by platform (long-read samples use minimap2)
    ch_reads
        .branch { meta, reads ->
            long_read:  meta.platform in ['ont', 'pacbio']
            short_read: true
        }
        .set { ch_branched }

    BWAMEM2_MEM (
        ch_branched.short_read,
        ch_reference.bwamem2_index,
        ch_reference.bwa_alt,
        ch_reference.fasta,
        true // coordinate-sort inside the aligner task
    )
    MINIMAP2_ALIGN (
        ch_branched.long_read,
        ch_reference.fasta
    )
    ch_sorted   = BWAMEM2_MEM.out.bam.mix(MINIMAP2_ALIGN.out.bam)
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first(), MINIMAP2_ALIGN.out.versions.first())

    SAMTOOLS_INDEX ( ch_sorted )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    bam      = ch_sorted.join(SAMTOOLS_INDEX.out.bai) // [ meta, bam, bai ]
    versions = ch_versions
}
