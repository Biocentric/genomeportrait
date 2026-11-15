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
include { DOWNLOAD_RESOURCE as DL_BWA_ALT } from '../../modules/local/download_resource'
include { DOWNLOAD_RESOURCE as DL_DBSNP   } from '../../modules/local/download_resource'
include { DOWNLOAD_RESOURCE as DL_KNOWN   } from '../../modules/local/download_resource'
include { DOWNLOAD_RESOURCE as DL_MILLS   } from '../../modules/local/download_resource'
include { DOWNLOAD_RESOURCE as DL_GNOMAD  } from '../../modules/local/download_resource'
include { DOWNLOAD_RESOURCE as DL_CLINVAR } from '../../modules/local/download_resource'
include { DOWNLOAD_VEP_CACHE             } from '../../modules/local/download_vep_cache'
include { PGS_PANEL_DOWNLOAD            } from '../../modules/local/pgs_panel_download'
include { DOWNLOAD_RESOURCE as DL_STR    } from '../../modules/local/download_resource'
include { DOWNLOAD_RESOURCE as DL_HIPSTR } from '../../modules/local/download_resource'
include { SAMTOOLS_FAIDX                } from '../../modules/local/samtools_faidx'
include { GATK4_CREATESEQUENCEDICTIONARY } from '../../modules/local/gatk4_dict'
include { BWAMEM2_INDEX                 } from '../../modules/local/bwamem2_index'
include { BUILD_REPORT_CONTAINER        } from '../../modules/local/build_report_container'
include { BUILD_HIPSTR_CONTAINER        } from '../../modules/local/build_hipstr_container'
include { BUILD_T1K_HLA                 } from '../../modules/local/build_t1k_hla'
include { BUILD_T1K_KIR                 } from '../../modules/local/build_t1k_kir'
include { DOWNLOAD_YLEAF_DATA           } from '../../modules/local/download_yleaf_data'
include { BUILD_YLEAF_CONTAINER        } from '../../modules/local/build_yleaf_container'
include { BUILD_CYRIUS_CONTAINER       } from '../../modules/local/build_cyrius_container'
include { BUILD_HMMIX_CONTAINER        } from '../../modules/local/build_hmmix_container'
include { HMMIX_RESOURCES              } from '../../modules/local/hmmix_resources'
include { BUILD_PAINTING_CONTAINER     } from '../../modules/local/build_painting_container'
include { PAINTING_REFERENCE           } from '../../modules/local/painting_reference'

workflow PREPARE_GENOME {

    take:
    genome          // val: genome key (e.g. GRCh38)
    needs_alignment // val: does any samplesheet row start from FASTQ?

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
    // The T1K KIR reference (IPD-KIR), only when KIR typing runs
    //
    ch_kir_index = Channel.value([])
    if (!params.skip_kir) {
        BUILD_T1K_KIR ()
        ch_kir_index = BUILD_T1K_KIR.out.index
        ch_versions  = ch_versions.mix(BUILD_T1K_KIR.out.versions)
    }

    //
    // Yleaf's Y-SNP data files (absent from the biocontainer), only when haplogroups run
    //
    ch_yleaf_data  = Channel.value([])
    ch_yleaf_ready = Channel.value('na')
    if (!params.skip_haplogroups) {
        DOWNLOAD_YLEAF_DATA ()
        ch_yleaf_data = DOWNLOAD_YLEAF_DATA.out.data
        ch_versions   = ch_versions.mix(DOWNLOAD_YLEAF_DATA.out.versions)
        if (params.build_yleaf_container) {
            BUILD_YLEAF_CONTAINER (
                file("${projectDir}/containers/yleaf/Singularity.def", checkIfExists: true)
            )
            ch_yleaf_ready = BUILD_YLEAF_CONTAINER.out.ready
            ch_versions    = ch_versions.mix(BUILD_YLEAF_CONTAINER.out.versions)
        }
    }

    //
    // Cyrius image (CYP2D6), only when pharmacogenomics runs with CYP2D6
    //
    ch_cyrius_ready = Channel.value('na')
    if (!params.skip_pharmcat && !params.skip_cyp2d6 && params.build_cyrius_container) {
        BUILD_CYRIUS_CONTAINER ( file("${projectDir}/containers/cyrius/Singularity.def", checkIfExists: true) )
        ch_cyrius_ready = BUILD_CYRIUS_CONTAINER.out.ready
        ch_versions     = ch_versions.mix(BUILD_CYRIUS_CONTAINER.out.versions)
    }

    //
    // hmmix image and resources (archaic segments), only when that step runs
    //
    ch_hmmix_ready     = Channel.value('na')
    ch_hmmix_resources = Channel.value([])
    if (!params.skip_archaic) {
        if (params.build_hmmix_container) {
            BUILD_HMMIX_CONTAINER ( file("${projectDir}/containers/hmmix/Singularity.def", checkIfExists: true) )
            ch_hmmix_ready = BUILD_HMMIX_CONTAINER.out.ready
            ch_versions    = ch_versions.mix(BUILD_HMMIX_CONTAINER.out.versions)
        }
        HMMIX_RESOURCES ( file("${projectDir}/scripts/zenodo_fetch.py", checkIfExists: true), ch_hmmix_ready )
        ch_hmmix_resources = HMMIX_RESOURCES.out.resources
        ch_versions        = ch_versions.mix(HMMIX_RESOURCES.out.versions)
    }

    //
    // Chromosome painting image and reference, only when painting runs
    //
    ch_painting_ready     = Channel.value('na')
    ch_painting_reference = Channel.value([])
    if (!params.skip_painting) {
        if (params.build_painting_container) {
            BUILD_PAINTING_CONTAINER ( file("${projectDir}/containers/painting/Singularity.def", checkIfExists: true) )
            ch_painting_ready = BUILD_PAINTING_CONTAINER.out.ready
            ch_versions       = ch_versions.mix(BUILD_PAINTING_CONTAINER.out.versions)
        }
        PAINTING_REFERENCE (
            file("${projectDir}/scripts/paint_reference.sh", checkIfExists: true),
            file("${projectDir}/assets/painting/holdout.tsv", checkIfExists: true),
            ch_painting_ready
        )
        ch_painting_reference = PAINTING_REFERENCE.out.reference
        ch_versions           = ch_versions.mix(PAINTING_REFERENCE.out.versions)
    }

    //
    // Core genome FASTA + indices
    //
    DL_FASTA ( [ [id:'genome'], refs.fasta ], "${params.reference_base}/${genome}/genome" )
    // A download dir can end up holding companion files (e.g. a stray .fai written by a tool
    // that resolved the symlink), and DOWNLOAD_RESOURCE emits a glob — so pick the FASTA
    // explicitly rather than passing whatever the glob matched.
    ch_fasta = DL_FASTA.out.file.map { meta, f ->
        (f instanceof List) ? f.find { it.name ==~ /(?i).*\.(fa|fasta|fna)(\.gz)?$/ } ?: f[0] : f
    }
    ch_versions = ch_versions.mix(DL_FASTA.out.versions)

    SAMTOOLS_FAIDX ( ch_fasta.map { [ [id:'genome'], it ] } )
    GATK4_CREATESEQUENCEDICTIONARY ( ch_fasta.map { [ [id:'genome'], it ] } )
    ch_fai  = SAMTOOLS_FAIDX.out.fai.map { it[1] }
    ch_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict.map { it[1] }
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    // bwa .alt file: ALT-aware mapping keeps reads at loci with ALT haplotypes (CYP2D6, HLA,
    // KIR, ...) on the primary assembly instead of scattering them onto ALT contigs at MAPQ 0.
    ch_bwa_alt = Channel.value([])
    if (needs_alignment && params.alt_aware && refs.bwa_alt) {
        DL_BWA_ALT ( [ [id:'bwa_alt'], refs.bwa_alt ], "${params.reference_base}/${genome}/bwa_alt" )
        ch_bwa_alt = DL_BWA_ALT.out.file.map { it[1] }
    }

    // bwa-mem2 index is large/slow; only build for short-read alignment from FASTQ
    ch_bwamem2 = Channel.value([])
    if (params.aligner == 'bwa-mem2' && needs_alignment) {
        BWAMEM2_INDEX ( ch_fasta.map { [ [id:'genome'], it ] } )
        ch_bwamem2  = BWAMEM2_INDEX.out.index
        ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)
    }

    //
    // Known sites for BQSR
    //
    // Only BQSR and HaplotypeCaller read these (dbSNP alone is ~11 GB); with the v2 defaults
    // (DeepVariant, no BQSR) nothing does, so they are not fetched.
    ch_dbsnp        = Channel.value([])
    ch_known_indels = Channel.value([])
    ch_mills        = Channel.value([])
    if (!params.skip_bqsr || params.variant_caller == 'haplotypecaller') {
        DL_DBSNP ( [ [id:'dbsnp'], refs.dbsnp ], "${params.reference_base}/${genome}/knownsites" )
        DL_KNOWN ( [ [id:'known_indels'], refs.known_indels ], "${params.reference_base}/${genome}/knownsites" )
        DL_MILLS ( [ [id:'mills'], refs.mills ], "${params.reference_base}/${genome}/knownsites" )
        ch_dbsnp        = DL_DBSNP.out.file.map { it[1] }
        ch_known_indels = DL_KNOWN.out.file.map { it[1] }
        ch_mills        = DL_MILLS.out.file.map { it[1] }
    }

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
    // One reference panel, HGDP + 1000 Genomes (the pgsc_calc release), for both ancestry and
    // polygenic scores: [ pgen, pvar.zst, psam, KING list of relatives to drop ].
    //
    ch_panel = Channel.value([])
    if (!params.skip_ancestry || !params.skip_prs) {
        if (params.reference_panel) {
            def d = params.reference_panel
            ch_panel = Channel.value([
                file("${d}/GRCh38_HGDP+1kGP_ALL.pgen",            checkIfExists: true),
                file("${d}/GRCh38_HGDP+1kGP_ALL.pvar.zst",        checkIfExists: true),
                file("${d}/GRCh38_HGDP+1kGP_ALL.psam",            checkIfExists: true),
                file("${d}/GRCh38_HGDP+1kGP.king.cutoff.out.id",  checkIfExists: true)
            ])
        } else {
            PGS_PANEL_DOWNLOAD ( params.reference_panel_url )
            ch_panel    = PGS_PANEL_DOWNLOAD.out.panel
            ch_versions = ch_versions.mix(PGS_PANEL_DOWNLOAD.out.versions)
        }
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
    bwa_alt       = ch_bwa_alt
    dbsnp         = ch_dbsnp
    known_indels  = ch_known_indels
    mills         = ch_mills
    vep_cache     = ch_vep_cache
    gnomad        = ch_gnomad
    clinvar       = ch_clinvar
    panel         = ch_panel
    str_catalog   = ch_str_catalog
    hipstr_codis  = ch_hipstr_codis
    report_sif    = ch_report_sif
    hipstr_ready  = ch_hipstr_ready
    hla_index     = ch_hla_index
    kir_index     = ch_kir_index
    yleaf_data    = ch_yleaf_data
    yleaf_ready   = ch_yleaf_ready
    cyrius_ready  = ch_cyrius_ready
    hmmix_ready   = ch_hmmix_ready
    hmmix_resources = ch_hmmix_resources
    painting_ready  = ch_painting_ready
    painting_reference = ch_painting_reference
    versions      = ch_versions
}
