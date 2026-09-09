//
// BAM_FORENSIC_STR: CODIS core STR profile (HipSTR) + pathogenic-locus repeat sizing (ExpansionHunter)
//   The CODIS profile is the same set of markers used in forensic identity databases.
//
include { CRAM_TO_BAM     } from '../../modules/local/cram_to_bam'
include { HIPSTR          } from '../../modules/local/hipstr'
include { EXPANSIONHUNTER } from '../../modules/local/expansionhunter'
include { MERGE_STR_CATALOG } from '../../modules/local/merge_str_catalog'
include { STR_REPORT      } from '../../modules/local/str_report'

workflow BAM_FORENSIC_STR {

    take:
    ch_bam       // channel: [ meta, cram, crai ]
    ch_reference // map

    main:
    ch_versions = Channel.empty()

    // HipSTR (and ExpansionHunter) are happiest with BAM; the analysis file is a CRAM,
    // and HipSTR's older htslib reads 0 reads from it. Convert once and feed both.
    CRAM_TO_BAM ( ch_bam, ch_reference.fasta, ch_reference.fai )
    ch_strbam = CRAM_TO_BAM.out.bam
    ch_versions = ch_versions.mix(CRAM_TO_BAM.out.versions.first())

    HIPSTR ( ch_strbam, ch_reference.fasta, ch_reference.fai, ch_reference.hipstr_codis, ch_reference.hipstr_ready )

    // HipSTR cannot reach CODIS loci whose reference tract exceeds what a 2x150 bp pair can
    // span (>71 bp: D21S11, D2S1338, FGA, vWA, D12S391) and drops them silently. Add them to
    // ExpansionHunter's catalog, which genotypes from a sequence graph using spanning,
    // flanking AND in-repeat reads, so it is not capped by read length.
    MERGE_STR_CATALOG ( ch_reference.str_catalog, file(params.codis_long_catalog, checkIfExists: true) )
    EXPANSIONHUNTER ( ch_strbam, ch_reference.fasta, ch_reference.fai, MERGE_STR_CATALOG.out.catalog )
    ch_versions = ch_versions.mix(HIPSTR.out.versions.first(), MERGE_STR_CATALOG.out.versions,
                                  EXPANSIONHUNTER.out.versions.first())

    STR_REPORT ( HIPSTR.out.vcf.join(EXPANSIONHUNTER.out.vcf) )
    ch_versions = ch_versions.mix(STR_REPORT.out.versions.first())

    emit:
    results  = STR_REPORT.out.tsv
    versions = ch_versions
}
