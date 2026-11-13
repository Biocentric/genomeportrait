//
// VCF_ANNOTATE: Ensembl VEP + SnpSift (gnomAD allele frequencies, ClinVar, dbNSFP)
//   Note: ClinVar is annotated for context only; the pipeline does not emit
//   disease-risk interpretations (see docs/ethics.md).
//
include { ENSEMBLVEP_VEP            } from '../../modules/local/ensemblvep_vep'
include { SNPSIFT_ANNOTATE as SNPSIFT_GNOMAD  } from '../../modules/local/snpsift_annotate'
include { SNPSIFT_ANNOTATE as SNPSIFT_CLINVAR } from '../../modules/local/snpsift_annotate'
include { VCF_TO_TABLE             } from '../../modules/local/vcf_to_table'

workflow VCF_ANNOTATE {

    take:
    ch_vcf       // channel: [ meta, vcf, tbi ]
    ch_reference // map

    main:
    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()

    ENSEMBLVEP_VEP (
        ch_vcf,
        ch_reference.fasta,
        ch_reference.vep_cache,
        params.vep_species,
        params.vep_cache_version
    )
    ch_versions = ch_versions.mix(ENSEMBLVEP_VEP.out.versions.first())
    ch_reports  = ch_reports.mix(ENSEMBLVEP_VEP.out.report)

    // gnomAD allele frequencies are pulled from the VEP cache (--af_gnomadg), so no separate
    // genome-wide gnomAD file is required. Optionally layer a gnomAD VCF in via SnpSift.
    ch_vep = ENSEMBLVEP_VEP.out.vcf
    if (params.use_snpsift_gnomad) {
        SNPSIFT_GNOMAD ( ch_vep, ch_reference.gnomad, 'gnomad' )
        ch_vep      = SNPSIFT_GNOMAD.out.vcf
        ch_versions = ch_versions.mix(SNPSIFT_GNOMAD.out.versions.first())
    }
    SNPSIFT_CLINVAR ( ch_vep, ch_reference.clinvar, 'clinvar' )

    // Flatten to a tidy TSV used by the report (one row per consequence)
    VCF_TO_TABLE ( SNPSIFT_CLINVAR.out.vcf )
    ch_versions = ch_versions.mix(VCF_TO_TABLE.out.versions.first())

    emit:
    vcf      = SNPSIFT_CLINVAR.out.vcf
    tsv      = VCF_TO_TABLE.out.tsv
    reports  = ch_reports
    versions = ch_versions
}
