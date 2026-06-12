//
// FASTQ_ALIGN: align short reads (bwa-mem2) or long reads (minimap2), sort + index
//
include { BWAMEM2_MEM    } from '../../modules/local/bwamem2_mem'
include { MINIMAP2_ALIGN } from '../../modules/local/minimap2_align'
include { SAMTOOLS_SORT  } from '../../modules/local/samtools_sort'
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
        ch_reference.fasta,
        true // sort
    )
    MINIMAP2_ALIGN (
        ch_branched.long_read,
        ch_reference.fasta
    )

    ch_bam_unsorted = BWAMEM2_MEM.out.bam.mix(MINIMAP2_ALIGN.out.bam)
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first(), MINIMAP2_ALIGN.out.versions.first())

    SAMTOOLS_SORT ( ch_bam_unsorted, ch_reference.fasta )
    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    ch_indexed = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai)

    emit:
    bam      = ch_indexed // [ meta, bam, bai ]
    versions = ch_versions
}
