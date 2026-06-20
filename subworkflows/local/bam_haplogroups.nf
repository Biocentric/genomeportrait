//
// BAM_HAPLOGROUPS: maternal (mtDNA, HaploGrep2) + paternal (Y, Yleaf) lineage
//
include { MITO_CALL    } from '../../modules/local/mito_call'
include { HAPLOGREP2   } from '../../modules/local/haplogrep2'
include { YLEAF        } from '../../modules/local/yleaf'
include { LINEAGE_REPORT } from '../../modules/local/lineage_report'

workflow BAM_HAPLOGROUPS {

    take:
    ch_bam       // channel: [ meta, cram, crai ]
    ch_vcf       // channel: [ meta, vcf, tbi ]
    ch_reference // map

    main:
    ch_versions = Channel.empty()

    // mtDNA haplogroup: call MT contig then classify with HaploGrep2 (PhyloTree 17)
    MITO_CALL ( ch_bam, ch_reference.fasta, ch_reference.fai )
    HAPLOGREP2 ( MITO_CALL.out.vcf )
    ch_versions = ch_versions.mix(MITO_CALL.out.versions.first(), HAPLOGREP2.out.versions.first())

    // Y haplogroup with Yleaf (only meaningful for XY samples; module no-ops otherwise)
    ch_bam_xy = ch_bam.filter { meta, bam, bai -> meta.sex != 'XX' }
    YLEAF ( ch_bam_xy, ch_reference.fasta, ch_reference.fai )
    ch_versions = ch_versions.mix(YLEAF.out.versions.first())

    // For XX samples Yleaf is skipped, so the remainder-join yields a null Y result;
    // convert it to an empty placeholder (null is rejected by the module's path input).
    LINEAGE_REPORT (
        HAPLOGREP2.out.result.join(YLEAF.out.hg, remainder: true)
            .map { meta, mt, y -> [ meta, mt, y ?: [] ] }
    )
    ch_versions = ch_versions.mix(LINEAGE_REPORT.out.versions.first())

    emit:
    results  = LINEAGE_REPORT.out.tsv
    versions = ch_versions
}
