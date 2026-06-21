/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREPARE_GENOME
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Download (once) and index every reference resource needed by the pipeline, and
    build the local report container. Everything is staged into
    params.reference_base/<genome>/ and re-used on subsequent runs (storeDir cache).

    Each resource is emitted as its own value channel, so the caller can
    assemble a plain map of value-channels (see main.nf). That keeps proper data
    dependencies while letting downstream code reference e.g. `ch_reference.fasta`.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DOWNLOAD_RESOURCE as DL_FASTA   } from '../../modules/local/download_resource'
include { DOWNLOAD_RESOURCE as DL_DBSNP   } from '../../modules/local/download_resource'
include { DOWNLOAD_RESOURCE as DL_KNOWN   } from '../../modules/local/download_resource'
include { DOWNLOAD_RESOURCE as DL_MILLS   } from '../../modules/local/download_resource'
include { DOWNLOAD_RESOURCE as DL_GNOMAD  } from '../../modules/local/download_resource'
include { DOWNLOAD_RESOURCE as DL_CLINVAR } from '../../modules/local/download_resource'
include { DOWNLOAD_VEP_CACHE             } from '../../modules/local/download_vep_cache'
include { DOWNLOAD_KGP_HGDP             } from '../../modules/local/download_kgp_hgdp'
include { DOWNLOAD_RESOURCE as DL_STR    } from '../../modules/local/download_resource'
include { DOWNLOAD_RESOURCE as DL_HIPSTR } from '../../modules/local/download_resource'
include { SAMTOOLS_FAIDX                } from '../../modules/local/samtools_faidx'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../modules/local/gatk4_dict'
include { BWAMEM2_INDEX                 } from '../../modules/local/bwamem2_index'
include { BUILD_REPORT_CONTAINER        } from '../../modules/local/build_report_container'
include { BUILD_HIPSTR_CONTAINER        } from '../../modules/local/build_hipstr_container'
include { BUILD_T1K_HLA                 } from '../../modules/local/build_t1k_hla'

workflow PREPARE_GENOME {

    take:
    genome // val: genome key (e.g. GRCh38)

    main:
    ch_versions = Channel.empty()
    // Accept common aliases for the human GRCh38/hg38 build
    if (genome.toLowerCase() in ['hg38', 'grch38', 'grch38_full', 'homo_sapiens'])
        genome = 'GRCh38'
    def refs = params.references[ genome ]
    if (refs == null)
        error "Unknown --genome '${genome}'. Known: ${params.references.keySet().join(', ')}"

    //
    // Build the local report container (cached). Other engines short-circuit inside.
    //
    ch_report_sif = Channel.empty()
    if (params.build_report_container) {
        BUILD_REPORT_CONTAINER (
            file("${projectDir}/containers/report/Singularity.def", checkIfExists: true),
            file("${projectDir}/containers/report/environment.yml", checkIfExists: true)
        )
        ch_report_sif = BUILD_REPORT_CONTAINER.out.ready
        ch_versions   = ch_versions.mix(BUILD_REPORT_CONTAINER.out.versions)
    }

    //
    // Build the local HipSTR container (compiled from source; only when the STR stage runs)
    //
    ch_hipstr_ready = Channel.value('na')
    if (!params.skip_str && params.build_hipstr_container) {
        BUILD_HIPSTR_CONTAINER (
            file("${projectDir}/containers/hipstr/Singularity.def", checkIfExists: true)
        )
        ch_hipstr_ready = BUILD_HIPSTR_CONTAINER.out.ready
        ch_versions     = ch_versions.mix(BUILD_HIPSTR_CONTAINER.out.versions)
    }

    //
    // Build the T1K HLA reference (IPD-IMGT/HLA), only when HLA typing runs
    //
    ch_hla_index = Channel.value([])
    if (!params.skip_hla) {
        BUILD_T1K_HLA ()
        ch_hla_index = BUILD_T1K_HLA.out.index
        ch_versions  = ch_versions.mix(BUILD_T1K_HLA.out.versions)
    }

    //
    // Core genome FASTA + indices
    //
    DL_FASTA ( [ [id:'genome'], refs.fasta ], "${params.reference_base}/${genome}/genome" )
    ch_fasta = DL_FASTA.out.file.map { meta, f -> f }
    ch_versions = ch_versions.mix(DL_FASTA.out.versions)

    SAMTOOLS_FAIDX ( ch_fasta.map { [ [id:'genome'], it ] } )
    GATK4_CREATESEQUENCEDICTIONARY ( ch_fasta.map { [ [id:'genome'], it ] } )
    ch_fai  = SAMTOOLS_FAIDX.out.fai.map { it[1] }
    ch_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict.map { it[1] }
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    // bwa-mem2 index is large/slow; only build for short-read alignment from FASTQ
    ch_bwamem2 = Channel.value([])
    if (params.aligner == 'bwa-mem2') {
        BWAMEM2_INDEX ( ch_fasta.map { [ [id:'genome'], it ] } )
        ch_bwamem2  = BWAMEM2_INDEX.out.index
        ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
    }

    //
    // Known sites for BQSR
    //
    DL_DBSNP ( [ [id:'dbsnp'], refs.dbsnp ], "${params.reference_base}/${genome}/knownsites" )
    DL_KNOWN ( [ [id:'known_indels'], refs.known_indels ], "${params.reference_base}/${genome}/knownsites" )
    DL_MILLS ( [ [id:'mills'], refs.mills ], "${params.reference_base}/${genome}/knownsites" )
    ch_dbsnp        = DL_DBSNP.out.file.map { it[1] }
    ch_known_indels = DL_KNOWN.out.file.map { it[1] }
    ch_mills        = DL_MILLS.out.file.map { it[1] }

    //
    // Annotation databases (only downloaded when annotation is enabled — the VEP cache is ~25 GB)
    //
    ch_vep_cache = Channel.value([])
    ch_gnomad    = Channel.value([])
    ch_clinvar   = Channel.value([])
    if (!params.skip_annotation) {
        DOWNLOAD_VEP_CACHE ( refs.vep_cache, "${params.reference_base}/${genome}/vep_cache" )
        DL_CLINVAR ( [ [id:'clinvar'], refs.clinvar ], "${params.reference_base}/${genome}/annotation" )
        ch_vep_cache = DOWNLOAD_VEP_CACHE.out.cache
        ch_clinvar   = DL_CLINVAR.out.file.map { it[1] }
        ch_versions  = ch_versions.mix(DOWNLOAD_VEP_CACHE.out.versions)
        // gnomAD VCF only needed if explicitly annotating with SnpSift (VEP cache covers AF otherwise)
        if (params.use_snpsift_gnomad) {
            DL_GNOMAD ( [ [id:'gnomad'], refs.gnomad ], "${params.reference_base}/${genome}/annotation" )
            ch_gnomad = DL_GNOMAD.out.file.map { it[1] }
        }
    }

    //
    // Ancestry panel (only if requested — it is large)
    //
    ch_kgp = Channel.value([])
    if (!params.skip_ancestry) {
        DOWNLOAD_KGP_HGDP (
            [ refs.kgp_hgdp_panel, refs.kgp_hgdp_pvar, refs.kgp_hgdp_psam ],
            "${params.reference_base}/${genome}/ancestry_panel"
        )
        ch_kgp = DOWNLOAD_KGP_HGDP.out.panel
        ch_versions = ch_versions.mix(DOWNLOAD_KGP_HGDP.out.versions)
    }

    //
    // STR catalogs (only when the forensic-STR stage is enabled)
    //
    ch_str_catalog  = Channel.value([])
    ch_hipstr_codis = Channel.value([])
    if (!params.skip_str) {
        DL_STR    ( [ [id:'str_catalog'],  refs.str_catalog  ], "${params.reference_base}/${genome}/str" )
        DL_HIPSTR ( [ [id:'hipstr_codis'], refs.hipstr_codis ], "${params.reference_base}/${genome}/str" )
        ch_str_catalog  = DL_STR.out.file.map { it[1] }
        ch_hipstr_codis = DL_HIPSTR.out.file.map { it[1] }
    }

    emit:
    fasta         = ch_fasta
    fai           = ch_fai
    dict          = ch_dict
    bwamem2_index = ch_bwamem2
    dbsnp         = ch_dbsnp
    known_indels  = ch_known_indels
    mills         = ch_mills
    vep_cache     = ch_vep_cache
    gnomad        = ch_gnomad
    clinvar       = ch_clinvar
    kgp_hgdp      = ch_kgp
    str_catalog   = ch_str_catalog
    hipstr_codis  = ch_hipstr_codis
    report_sif    = ch_report_sif
    hipstr_ready  = ch_hipstr_ready
    hla_index     = ch_hla_index
    versions      = ch_versions
}
