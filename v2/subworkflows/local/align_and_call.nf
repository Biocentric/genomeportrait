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
include { PARABRICKS_FQ2BAM_DOCKER } from '../../modules/local/parabricks_fq2bam_docker'
include { PARABRICKS_DEEPVARIANT   } from '../../modules/nf-core/parabricks/deepvariant/main'
include { BCFTOOLS_NORM as BCFTOOLS_NORM_GPU } from '../../modules/local/bcftools_norm'

include { BWAMEM2_MARKDUP_CRAM     } from '../../modules/local/bwamem2_markdup_cram'
include { FASTQ_ALIGN              } from './fastq_align'
include { BAM_MARKDUP_BQSR         } from './bam_markdup_bqsr'
include { BAM_CALL_SMALLVARIANTS   } from './bam_call_smallvariants'

workflow ALIGN_AND_CALL {

    take:
    ch_reads       // channel: [ meta, [ reads ] ]   (trimmed FASTQ)
    ch_prealigned  // channel: [ meta, bam/cram, index ] (analysis-ready samplesheet alignments without calls)
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

    // The BQSR recalibration table fq2bam can emit is only consumed when the GPU
    // caller runs. Otherwise skip --knownSites: it would stage a ~11 GB dbSNP VCF
    // into the task dir (stageInMode 'copy') to produce a table nothing reads.
    ch_gpu_known_sites = params.gpu_variant_calling && !params.skip_bqsr
        ? ch_reference.dbsnp.map { [ [id:'dbsnp'], it ] }
        : Channel.value([ [id:'dbsnp'], [] ])

    ch_gpu_reads = ch_routed.gpu.map { meta, reads -> [ meta, reads, [] ] }
    ch_gpu_fasta = ch_reference.fasta.map { [ [id:'genome'], it ] }

    if (params.parabricks_launcher == 'docker') {
        PARABRICKS_FQ2BAM_DOCKER ( ch_gpu_reads, ch_gpu_fasta, BWA_INDEX.out.index, ch_gpu_known_sites )
        ch_gpu_bam  = PARABRICKS_FQ2BAM_DOCKER.out.bam.join(PARABRICKS_FQ2BAM_DOCKER.out.bai)
        ch_versions = ch_versions.mix(PARABRICKS_FQ2BAM_DOCKER.out.versions)
    }
    else {
        PARABRICKS_FQ2BAM ( ch_gpu_reads, ch_gpu_fasta, BWA_INDEX.out.index, ch_gpu_known_sites )
        ch_gpu_bam  = PARABRICKS_FQ2BAM.out.bam.join(PARABRICKS_FQ2BAM.out.bai)
        ch_versions = ch_versions.mix(PARABRICKS_FQ2BAM.out.versions)
    }

    // pbrun deepvariant requires compute capability >= 7.0 even where fq2bam does not.
    // On Pascal (P40/P100 + Parabricks 4.0) set --gpu_variant_calling false: alignment
    // still runs on the GPU and the Parabricks BAM is handed to the CPU caller below,
    // which also restores the gVCF that the GPU path cannot produce.
    ch_gpu_vcf     = Channel.empty()
    ch_gpu_for_cpu = Channel.empty()
    if (params.gpu_variant_calling) {
        PARABRICKS_DEEPVARIANT (
            ch_gpu_bam.map { meta, bam, bai -> [ meta, bam, bai, [] ] },
            ch_reference.fasta.map { [ [id:'genome'], it ] }
        )
        ch_versions = ch_versions.mix(PARABRICKS_DEEPVARIANT.out.versions)

        // Normalise the GPU VCF to the same shape as the CPU caller output [meta, vcf.gz, tbi]
        BCFTOOLS_NORM_GPU ( PARABRICKS_DEEPVARIANT.out.vcf.map { m, v -> [ m, v, [] ] }, ch_reference.fasta )
        ch_gpu_vcf  = BCFTOOLS_NORM_GPU.out.vcf
        ch_versions = ch_versions.mix(BCFTOOLS_NORM_GPU.out.versions)
    }
    else {
        ch_gpu_for_cpu = ch_gpu_bam
    }

    //
    // ---------------- CPU path (default toolchain) ----------------
    //
    // Short reads with bwa-mem2 go straight to a duplicate-marked CRAM in one stream
    // (--markdup samtools, the default). Long reads, minimap2, and --markdup gatk keep the
    // separate path: sorted BAM, then GATK MarkDuplicates.
    ch_routed.cpu
        .branch { meta, reads ->
            stream: params.markdup == 'samtools' && params.aligner == 'bwa-mem2' && !(meta.platform in ['ont', 'pacbio'])
            classic: true
        }
        .set { ch_cpu }
    BWAMEM2_MARKDUP_CRAM ( ch_cpu.stream, ch_reference.bwamem2_index, ch_reference.bwa_alt, ch_reference.fasta, ch_reference.fai )
    FASTQ_ALIGN ( ch_cpu.classic, ch_reference )
    BAM_MARKDUP_BQSR ( FASTQ_ALIGN.out.bam, ch_reference, params.skip_bqsr, BWAMEM2_MARKDUP_CRAM.out.cram )
    // A samplesheet alignment is taken as analysis-ready (sorted, duplicates marked): it is
    // called as it is. v1 pushed it through MarkDuplicates and BQSR a second time.
    BAM_CALL_SMALLVARIANTS (
        BAM_MARKDUP_BQSR.out.bam.mix(ch_gpu_for_cpu).mix(ch_prealigned),
        ch_reference
    )
    ch_versions = ch_versions.mix(
        BWAMEM2_MARKDUP_CRAM.out.versions.first(),
        FASTQ_ALIGN.out.versions,
        BAM_MARKDUP_BQSR.out.versions,
        BAM_CALL_SMALLVARIANTS.out.versions
    )

    //
    // ---------------- Merge GPU + CPU outputs ----------------
    //
    emit:
    bam      = BAM_MARKDUP_BQSR.out.bam.mix(ch_gpu_bam)        // [ meta, bam/cram, bai/crai ]
    // A real BAM for tools that cannot read CRAM (Qualimap): the MarkDuplicates BAM on the
    // separate path, the fq2bam BAM on the GPU path. The streamed path writes no BAM.
    md_bam   = BAM_MARKDUP_BQSR.out.md_bam.mix(ch_gpu_bam)     // [ meta, bam, bai ]
    vcf      = BAM_CALL_SMALLVARIANTS.out.vcf.mix(ch_gpu_vcf)  // [ meta, vcf.gz, tbi ]
    // gVCF (hom-ref blocks) is required by VCF_POLYGENIC. CPU path only: `pbrun
    // deepvariant` is not invoked with --gvcf, so GPU samples emit nothing here and
    // the polygenic stage simply does not run for them.
    gvcf     = BAM_CALL_SMALLVARIANTS.out.gvcf                 // [ meta, gvcf.gz, tbi ]
    metrics  = BAM_MARKDUP_BQSR.out.metrics.mix(BWAMEM2_MARKDUP_CRAM.out.metrics)
    stats    = BAM_CALL_SMALLVARIANTS.out.stats
    versions = ch_versions
}
