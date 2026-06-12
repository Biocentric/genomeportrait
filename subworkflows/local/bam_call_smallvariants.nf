//
// BAM_CALL_SMALLVARIANTS: germline SNV/indel calling
//   default = DeepVariant; alternatives = GATK HaplotypeCaller, Clair3 (long read)
//
include { DEEPVARIANT_RUNDEEPVARIANT } from '../../modules/local/deepvariant'
include { GATK4_HAPLOTYPECALLER      } from '../../modules/local/gatk4_haplotypecaller'
include { CLAIR3                     } from '../../modules/local/clair3'
include { BCFTOOLS_NORM             } from '../../modules/local/bcftools_norm'
include { BCFTOOLS_STATS            } from '../../modules/local/bcftools_stats'

workflow BAM_CALL_SMALLVARIANTS {

    take:
    ch_bam       // channel: [ meta, cram, crai ]
    ch_reference // map

    main:
    ch_versions = Channel.empty()
    ch_raw_vcf  = Channel.empty()
    ch_gvcf     = Channel.empty()

    if (params.variant_caller == 'deepvariant') {
        DEEPVARIANT_RUNDEEPVARIANT ( ch_bam, ch_reference.fasta, ch_reference.fai )
        ch_raw_vcf  = DEEPVARIANT_RUNDEEPVARIANT.out.vcf
        ch_gvcf     = DEEPVARIANT_RUNDEEPVARIANT.out.gvcf
        ch_versions = ch_versions.mix(DEEPVARIANT_RUNDEEPVARIANT.out.versions.first())
    } else if (params.variant_caller == 'clair3') {
        CLAIR3 ( ch_bam, ch_reference.fasta, ch_reference.fai )
        ch_raw_vcf  = CLAIR3.out.vcf
        ch_versions = ch_versions.mix(CLAIR3.out.versions.first())
    } else {
        GATK4_HAPLOTYPECALLER ( ch_bam, ch_reference.fasta, ch_reference.fai, ch_reference.dict, ch_reference.dbsnp )
        ch_raw_vcf  = GATK4_HAPLOTYPECALLER.out.vcf
        ch_gvcf     = GATK4_HAPLOTYPECALLER.out.gvcf
        ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions.first())
    }

    // Left-align + split multiallelics for downstream consistency
    BCFTOOLS_NORM ( ch_raw_vcf, ch_reference.fasta )
    BCFTOOLS_STATS ( BCFTOOLS_NORM.out.vcf )
    ch_versions = ch_versions.mix(BCFTOOLS_NORM.out.versions.first())

    emit:
    vcf      = BCFTOOLS_NORM.out.vcf    // [ meta, vcf.gz, tbi ]
    gvcf     = ch_gvcf
    stats    = BCFTOOLS_STATS.out.stats
    versions = ch_versions
}
