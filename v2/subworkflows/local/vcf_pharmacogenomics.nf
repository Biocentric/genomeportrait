//
// BAM_PHARMACOGENOMICS: PharmCAT star-allele calling + drug-gene annotations
//   Reports metaboliser phenotypes (e.g. CYP2D6, CYP2C19) and CPIC guidance.
//   Informational only — not a prescription tool.
//
include { CYRIUS          } from '../../modules/local/cyrius'
include { PHARMCAT        } from '../../modules/local/pharmcat'
include { PHARMCAT_REPORT } from '../../modules/local/pharmcat_report'

workflow BAM_PHARMACOGENOMICS {

    take:
    ch_bam       // channel: [ meta, bam/cram, bai/crai ]
    ch_reference // map

    main:
    ch_versions = Channel.empty()

    // CYP2D6 from the alignment (Cyrius), handed to PharmCAT as an outside call
    if (!params.skip_cyp2d6) {
        CYRIUS ( ch_bam, ch_reference.fasta, ch_reference.fai, ch_reference.cyrius_ready )
        ch_versions = ch_versions.mix(CYRIUS.out.versions.first())
        ch_pharmcat_in = ch_bam
            .map { meta, bam, bai -> [ meta.id, meta, bam, bai ] }
            .join( CYRIUS.out.outside.map { meta, f -> [ meta.id, f ] } )
            .map { id, meta, bam, bai, outside -> [ meta, bam, bai, outside ] }
    } else {
        ch_pharmcat_in = ch_bam.map { meta, bam, bai -> [ meta, bam, bai, [] ] }
    }
    PHARMCAT ( ch_pharmcat_in, ch_reference.fasta, ch_reference.fai )
    ch_versions = ch_versions.mix(PHARMCAT.out.versions.first())

    PHARMCAT_REPORT ( PHARMCAT.out.json, file("${projectDir}/scripts/pharmcat_report.py", checkIfExists: true) )
    ch_versions = ch_versions.mix(PHARMCAT_REPORT.out.versions.first())

    emit:
    results  = PHARMCAT_REPORT.out.tsv
    cyp2d6   = params.skip_cyp2d6 ? Channel.empty() : CYRIUS.out.tsv
    versions = ch_versions
}
