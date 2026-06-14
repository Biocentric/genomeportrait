//
// ALIGN_AND_CALL: align + mark-dup/BQSR + germline small-variant calling, routed at
// runtime to NVIDIA Parabricks (GPU) or the CPU toolchain based on CHECK_GPU.
//
//   GPU path (Parabricks): fq2bam (bwa-mem GPU + markdup [+ BQSR table]) -> DeepVariant GPU
//   CPU path (default):    bwa-mem2/minimap2 -> GATK markdup+BQSR -> DeepVariant/HC/Clair3
//
// Only one path receives data per sample (the other's processes simply don't execute),
// so the GPU container/index are never built on CPU-only runs.
//
include { BWA_INDEX                } from '../../modules/local/bwa_index'
include { PARABRICKS_FQ2BAM        } from '../../modules/nf-core/parabricks/fq2bam/main'
include { PARABRICKS_DEEPVARIANT   } from '../../modules/nf-core/parabricks/deepvariant/main'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_GPU } from '../../modules/local/bcftools_norm'

include { FASTQ_ALIGN              } from './fastq_align'
include { BAM_MARKDUP_BQSR         } from './bam_markdup_bqsr'
include { BAM_CALL_SMALLVARIANTS   } from './bam_call_smallvariants'

workflow ALIGN_AND_CALL {

    take:
    ch_reads       // channel: [ meta, [ reads ] ]   (trimmed FASTQ)
    ch_prealigned  // channel: [ meta, bam, bai ]     (samplesheet BAM/CRAM inputs)
    ch_reference   // map
    ch_gpu_ok      // value channel: Boolean — is a Parabricks-capable GPU available?

    main:
    ch_versions = Channel.empty()

    // Route FASTQ reads. Parabricks handles short reads only; long reads always go CPU.
    ch_reads.combine(ch_gpu_ok)
        .branch { meta, reads, ok ->
            gpu: ok && !(meta.platform in ['ont', 'pacbio'])
                return [ meta, reads ]
            cpu: true
                return [ meta, reads ]
        }
        .set { ch_routed }

    //
    // ---------------- GPU path (Parabricks) ----------------
    //   Build the classic BWA index only when the GPU path actually has reads.
    //
    ch_bwa_in = ch_routed.gpu
        .combine(ch_reference.fasta)
        .map { meta, reads, fasta -> [ [id:'genome'], fasta ] }
        .first()
    BWA_INDEX ( ch_bwa_in )
    ch_versions = ch_versions.mix(BWA_INDEX.out.versions)

    PARABRICKS_FQ2BAM (
        ch_routed.gpu.map { meta, reads -> [ meta, reads, [] ] },
        ch_reference.fasta.map { [ [id:'genome'], it ] },
        BWA_INDEX.out.index,
        ch_reference.dbsnp.map { [ [id:'dbsnp'], it ] }
    )
    ch_gpu_bam = PARABRICKS_FQ2BAM.out.bam.join(PARABRICKS_FQ2BAM.out.bai)
    ch_versions = ch_versions.mix(PARABRICKS_FQ2BAM.out.versions)

    PARABRICKS_DEEPVARIANT (
        ch_gpu_bam.map { meta, bam, bai -> [ meta, bam, bai, [] ] },
        ch_reference.fasta.map { [ [id:'genome'], it ] }
    )
    ch_versions = ch_versions.mix(PARABRICKS_DEEPVARIANT.out.versions)

    // Normalise the GPU VCF to the same shape as the CPU caller output [meta, vcf.gz, tbi]
    BCFTOOLS_NORM_GPU ( PARABRICKS_DEEPVARIANT.out.vcf.map { m, v -> [ m, v, [] ] }, ch_reference.fasta )
    ch_gpu_vcf = BCFTOOLS_NORM_GPU.out.vcf
    ch_versions = ch_versions.mix(BCFTOOLS_NORM_GPU.out.versions)

    //
    // ---------------- CPU path (default toolchain) ----------------
    //
    FASTQ_ALIGN ( ch_routed.cpu, ch_reference )
    ch_cpu_aligned = FASTQ_ALIGN.out.bam.mix(ch_prealigned)
    BAM_MARKDUP_BQSR ( ch_cpu_aligned, ch_reference, params.skip_bqsr )
    BAM_CALL_SMALLVARIANTS ( BAM_MARKDUP_BQSR.out.bam, ch_reference )
    ch_versions = ch_versions.mix(
        FASTQ_ALIGN.out.versions,
        BAM_MARKDUP_BQSR.out.versions,
        BAM_CALL_SMALLVARIANTS.out.versions
    )

    //
    // ---------------- Merge GPU + CPU outputs ----------------
    //
    emit:
    bam      = BAM_MARKDUP_BQSR.out.bam.mix(ch_gpu_bam)        // [ meta, bam/cram, bai/crai ]
    vcf      = BAM_CALL_SMALLVARIANTS.out.vcf.mix(ch_gpu_vcf)  // [ meta, vcf.gz, tbi ]
    metrics  = BAM_MARKDUP_BQSR.out.metrics
    stats    = BAM_CALL_SMALLVARIANTS.out.stats
    versions = ch_versions
}
