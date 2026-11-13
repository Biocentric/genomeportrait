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
    ch_sex       // channel: [ meta.id, 'XX' | 'XY' | 'unknown' ] — samplesheet value, else inferred from depth

    main:
    ch_versions = Channel.empty()

    // mtDNA haplogroup: call MT contig then classify with HaploGrep2 (PhyloTree 17)
    MITO_CALL ( ch_bam, ch_reference.fasta, ch_reference.fai )
    HAPLOGREP2 ( MITO_CALL.out.vcf )
    ch_versions = ch_versions.mix(MITO_CALL.out.versions.first(), HAPLOGREP2.out.versions.first())

    // Y haplogroup with Yleaf (only meaningful for XY samples; module no-ops otherwise)
    // v1 read meta.sex, so an undeclared sex always ran Yleaf — on XX samples too. The
    // resolved sex (GENOME_FACTS) decides; 'unknown' still runs, since Yleaf itself answers it.
    ch_bam_xy = ch_bam
        .map { meta, bam, bai -> [ meta.id, meta, bam, bai ] }
        .join(ch_sex)
        .filter { id, meta, bam, bai, sex -> sex != 'XX' }
        .map { id, meta, bam, bai, sex -> [ meta, bam, bai ] }
    YLEAF ( ch_bam_xy, ch_reference.fasta, ch_reference.fai, ch_reference.yleaf_data, ch_reference.yleaf_ready )
    ch_versions = ch_versions.mix(YLEAF.out.versions.first())

    // For XX samples Yleaf is skipped, so the remainder-join yields a null Y result;
    // convert it to an empty placeholder (null is rejected by the module's path input).
    LINEAGE_REPORT (
        HAPLOGREP2.out.result.join(YLEAF.out.hg, remainder: true)
            .map { meta, mt, y -> [ meta, mt, y ?: [] ] },
        file("${projectDir}/scripts/lineage_report.py", checkIfExists: true)
    )
    ch_versions = ch_versions.mix(LINEAGE_REPORT.out.versions.first())

    emit:
    results  = LINEAGE_REPORT.out.tsv
    versions = ch_versions
}
