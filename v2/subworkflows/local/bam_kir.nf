//
// BAM_KIR: KIR gene content and alleles with T1K (NK-cell receptor genes, chr19 LRC)
//
include { KIR_EXTRACT_FASTQ } from '../../modules/local/kir_extract_fastq'
include { KIR_T1K           } from '../../modules/local/kir_t1k'
include { KIR_REPORT        } from '../../modules/local/kir_report'

workflow BAM_KIR {

    take:
    ch_bam       // channel: [ meta, cram, crai ]
    ch_reference // map; .kir_index from PREPARE_GENOME

    main:
    ch_versions = Channel.empty()

    KIR_EXTRACT_FASTQ ( ch_bam, ch_reference.fasta, ch_reference.fai )
    KIR_T1K ( KIR_EXTRACT_FASTQ.out.reads, ch_reference.kir_index )
    KIR_REPORT ( KIR_T1K.out.alleles, file("${projectDir}/scripts/kir_report.py", checkIfExists: true) )
    ch_versions = ch_versions.mix(KIR_EXTRACT_FASTQ.out.versions.first(), KIR_T1K.out.versions.first(),
                                  KIR_REPORT.out.versions.first())

    emit:
    results  = KIR_REPORT.out.tsv
    versions = ch_versions
}
