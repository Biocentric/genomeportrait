//
// VCF_POLYGENIC: ancestry-adjusted polygenic scores with interpretation.
//
//   Sample-independent, cached under reference_base/pgs/<panel_key>/ (first genome pays):
//     PGS_FETCH            harmonised GRCh38 scoring files + catalogue metadata
//     PGS_SITES            union of scored positions
//     PGS_REFERENCE        HGDP+1kGP subset to those positions, unrelated only; PCA
//     PGS_MATCH            every score matched to the reference by chr:pos:ref:alt
//     PGS_SCORE_REFERENCE  the reference distribution of every score
//   Per sample (minutes):
//     PGS_TARGET           gVCF -> the sample in the reference's allele space
//     PGS_PROJECT          projection onto the reference PCs (plink2)
//     PGS_REPORT           single-pass scoring, ancestry-adjusted Z, direction, consensus, effect size, overview page
//
//   The quantities are those of pgsc_calc v2 (Z_norm2, percentile in the most similar
//   population); the interpretation is prs-selfscore's. See docs/polygenic.md.
//
include { PGS_FETCH           } from '../../modules/local/pgs_fetch'
include { PGS_PANEL_DOWNLOAD  } from '../../modules/local/pgs_panel_download'
include { PGS_SITES           } from '../../modules/local/pgs_sites'
include { PGS_REFERENCE       } from '../../modules/local/pgs_reference'
include { PGS_MATCH           } from '../../modules/local/pgs_match'
include { PGS_SCORE_REFERENCE } from '../../modules/local/pgs_score_reference'
include { PGS_TARGET          } from '../../modules/local/pgs_target'
include { PGS_PROJECT         } from '../../modules/local/pgs_project'
include { PGS_REPORT          } from '../../modules/local/pgs_report'

workflow VCF_POLYGENIC {

    take:
    ch_gvcf // channel: [ meta, gvcf, tbi ] — hom-ref blocks are required, a variant-only VCF will not do

    main:
    ch_versions = Channel.empty()

    def scripts   = "${projectDir}/scripts"
    def engine    = file("${scripts}/pgs_engine.py",    checkIfExists: true)
    def panel     = file(params.pgs_panel,              checkIfExists: true)
    def extra_ids = params.pgs_ids ? params.pgs_ids.toString().replaceAll(/\s/, '') : ''

    // Cache key: panel name + a digest of everything that changes the sample-independent
    // results (the score list and the ambiguous-SNP policy). Editing a concept or a note
    // in the panel table does not invalidate the cache; adding a score does.
    def ids = panel.readLines()
        .findAll { it && !it.startsWith('#') }
        .drop(1)
        .collect { it.tokenize('\t')[0] } + (extra_ids ? extra_ids.tokenize(',') : [])
    def panel_key = panel.baseName + '_' +
        (ids.unique().sort().join(',') + "|amb=${params.pgs_keep_ambiguous}|pcs=${params.pgs_n_pcs}").md5().take(10)

    PGS_FETCH (
        panel_key, panel, extra_ids,
        file(params.pgs_catalog, checkIfExists: true),
        file("${scripts}/pgs_fetch.py",   checkIfExists: true),
        file("${scripts}/pgs_catalog.py", checkIfExists: true)
    )
    ch_versions = ch_versions.mix(PGS_FETCH.out.versions)

    // Reference panel: a local, already-extracted copy if given, else download once.
    if (params.pgs_reference_panel) {
        def d = params.pgs_reference_panel
        ch_panel = Channel.value([
            file("${d}/GRCh38_HGDP+1kGP_ALL.pgen",            checkIfExists: true),
            file("${d}/GRCh38_HGDP+1kGP_ALL.pvar.zst",        checkIfExists: true),
            file("${d}/GRCh38_HGDP+1kGP_ALL.psam",            checkIfExists: true),
            file("${d}/GRCh38_HGDP+1kGP.king.cutoff.out.id",  checkIfExists: true)
        ])
    } else {
        PGS_PANEL_DOWNLOAD ( params.pgs_reference_panel_url )
        ch_panel    = PGS_PANEL_DOWNLOAD.out.panel
        ch_versions = ch_versions.mix(PGS_PANEL_DOWNLOAD.out.versions)
    }

    PGS_SITES ( panel_key, PGS_FETCH.out.scorefiles, engine )
    PGS_REFERENCE (
        panel_key, ch_panel, PGS_SITES.out.sites,
        file("${projectDir}/assets/pgs/pca_exclude_grch38.bed1", checkIfExists: true)
    )
    PGS_MATCH ( panel_key, PGS_FETCH.out.scorefiles, PGS_REFERENCE.out.pfile, PGS_REFERENCE.out.afreq, engine )
    PGS_SCORE_REFERENCE ( panel_key, PGS_REFERENCE.out.pfile, PGS_MATCH.out.scores )
    ch_versions = ch_versions.mix(
        PGS_SITES.out.versions, PGS_REFERENCE.out.versions,
        PGS_MATCH.out.versions, PGS_SCORE_REFERENCE.out.versions
    )

    PGS_TARGET ( ch_gvcf, PGS_REFERENCE.out.pfile, engine )
    PGS_PROJECT ( PGS_TARGET.out.vcf, PGS_REFERENCE.out.pca )
    PGS_REPORT (
        PGS_TARGET.out.vcf.join(PGS_PROJECT.out.pcs),
        PGS_SCORE_REFERENCE.out.scores,
        PGS_REFERENCE.out.ref_pcs,
        PGS_REFERENCE.out.pfile,
        PGS_REFERENCE.out.afreq,
        PGS_MATCH.out.scores,
        PGS_MATCH.out.stats,
        panel,
        file(params.pgs_concepts,      checkIfExists: true),
        file(params.pgs_r2_exclusions, checkIfExists: true),
        PGS_FETCH.out.catalog,
        engine,
        file("${scripts}/pgs_interpret.py", checkIfExists: true),
        file("${scripts}/pgs_overview.py",  checkIfExists: true)
    )
    ch_versions = ch_versions.mix(
        PGS_TARGET.out.versions.first(), PGS_PROJECT.out.versions.first(), PGS_REPORT.out.versions.first()
    )

    emit:
    fragment = PGS_REPORT.out.fragment   // [ meta, html fragment ] for the unified report
    tsv      = PGS_REPORT.out.tsv        // [ meta, interpreted tsv ]
    html     = PGS_REPORT.out.html
    versions = ch_versions
}
