/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// QC / prep
include { FASTQ_QC_TRIM            } from '../subworkflows/local/fastq_qc_trim'
include { CHECK_GPU                } from '../modules/local/check_gpu'
include { ALIGN_AND_CALL           } from '../subworkflows/local/align_and_call'
include { BAM_QC                   } from '../subworkflows/local/bam_qc'

// Variants
include { BAM_CALL_SV_CNV          } from '../subworkflows/local/bam_call_sv_cnv'
include { VCF_ANNOTATE             } from '../subworkflows/local/vcf_annotate'

// Interpretation layers
include { BAM_ANCESTRY             } from '../subworkflows/local/vcf_ancestry'
include { BAM_HAPLOGROUPS          } from '../subworkflows/local/bam_haplogroups'
include { BAM_FORENSIC_STR         } from '../subworkflows/local/bam_forensic_str'
include { VCF_POLYGENIC            } from '../subworkflows/local/vcf_polygenic'
include { MOSDEPTH                 } from '../modules/local/mosdepth'
include { GENOME_FACTS             } from '../modules/local/genome_facts'
include { VERIFYBAMID2             } from '../modules/local/verifybamid2'
include { CONTAMINATION_REPORT     } from '../modules/local/contamination_report'
include { ARRAY_PANEL              } from '../modules/local/array_panel'
include { EXPORT_23ANDME           } from '../modules/local/export_23andme'
include { RHD_DEPTH                } from '../modules/local/rhd_depth'
include { BLOOD_GROUPS             } from '../modules/local/blood_groups'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_INPUT } from '../modules/local/bcftools_stats'
include { BAM_HLA                  } from '../subworkflows/local/bam_hla'
include { BAM_KIR                  } from '../subworkflows/local/bam_kir'
include { ARCHAIC_HMMIX            } from '../modules/local/archaic_hmmix'
include { PAINT_CHROMOSOMES        } from '../modules/local/paint_chromosomes'
include { BAM_PHARMACOGENOMICS     } from '../subworkflows/local/vcf_pharmacogenomics'
include { VCF_TRAITS               } from '../subworkflows/local/vcf_traits'
include { BAM_EXTRAS               } from '../subworkflows/local/bam_extras'

// Reporting
include { GENOMEPORTRAIT_REPORT    } from '../modules/local/genomeportrait_report'
include { MULTIQC                  } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap         } from 'plugin/nf-schema'
include { paramsSummaryMultiqc     } from '../subworkflows/local/utils_genomeportrait_pipeline'
include { softwareVersionsToYAML   } from '../subworkflows/local/utils_genomeportrait_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELPERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// A samplesheet cell that was left empty arrives as null or as an empty list.
def isGiven(x) {
    return x != null && x != []
}

// Index next to a samplesheet file: <f>.crai / .bai / .tbi, or the extension swapped.
def withIndex(f) {
    def n = f.toString()
    def cands = n.endsWith('.cram') ? [ n + '.crai', n.replaceAll(/\.cram$/, '.crai') ]
              : n.endsWith('.bam')  ? [ n + '.bai',  n.replaceAll(/\.bam$/,  '.bai')  ]
              :                       [ n + '.tbi',  n + '.csi' ]
    def idx = cands.collect { c -> file(c) }.find { c -> c.exists() }
    if (!idx)
        error "No index found for ${n} (looked for ${cands.join(', ')})"
    return idx
}

// Hold a channel back until `ready` has emitted (the report container build).
def gateOn(ch, ready) {
    return ch.combine(ready).map { row -> row[0..-2] }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GENOMEPORTRAIT {

    take:
    ch_samplesheet // channel: [ meta, fastq_1, fastq_2, bam, gvcf, vcf ]
    ch_reference   // map    : prepared reference resources from PREPARE_GENOME

    main:

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()
    // Per-sample result collectors feeding the final portrait report
    ch_report_parts  = Channel.empty()

    //
    // Three ways in, decided per samplesheet row:
    //   fastq     reads only                    -> trim, align, mark duplicates, call, interpret
    //   aligned   bam/cram only                 -> call, interpret   (alignment taken as analysis-ready)
    //   portrait  bam/cram + gvcf + vcf         -> interpret only    (re-portray a finished genome)
    //
    ch_samplesheet
        .map { meta, fq1, fq2, bam, gvcf, vcf ->
            if (isGiven(gvcf) != isGiven(vcf))
                error "Sample ${meta.id}: give both 'gvcf' and 'vcf' (with 'bam') to skip variant calling, or neither."
            if ((isGiven(gvcf) || isGiven(vcf)) && !isGiven(bam))
                error "Sample ${meta.id}: 'gvcf'/'vcf' need the alignment in 'bam' too (STR, HLA, lineage and PGx read it)."
            def mode = isGiven(bam) ? (isGiven(gvcf) ? 'portrait' : 'aligned') : 'fastq'
            [ meta, mode, [ fq1, fq2 ].findAll { f -> isGiven(f) }, bam, gvcf, vcf ]
        }
        .branch { meta, mode, reads, bam, gvcf, vcf ->
            fastq:    mode == 'fastq'
            aligned:  mode == 'aligned'
            portrait: mode == 'portrait'
        }
        .set { ch_inputs }

    ch_prealigned    = ch_inputs.aligned.map  { meta, mode, reads, bam, gvcf, vcf -> [ meta, bam,  withIndex(bam)  ] }
    ch_portrait_bam  = ch_inputs.portrait.map { meta, mode, reads, bam, gvcf, vcf -> [ meta, bam,  withIndex(bam)  ] }
    ch_portrait_gvcf = ch_inputs.portrait.map { meta, mode, reads, bam, gvcf, vcf -> [ meta, gvcf, withIndex(gvcf) ] }
    ch_portrait_vcf  = ch_inputs.portrait.map { meta, mode, reads, bam, gvcf, vcf -> [ meta, vcf,  withIndex(vcf)  ] }

    if (!params.skip_bqsr && params.variant_caller == 'deepvariant')
        log.warn "[genomeportrait] BQSR is enabled with DeepVariant. DeepVariant is trained on unrecalibrated " +
                 "base qualities; this costs a whole-genome pass and does not improve its calls."

    //
    // STAGE 1: Read QC + trimming (FASTQ rows only)
    //
    FASTQ_QC_TRIM ( ch_inputs.fastq.map { meta, mode, reads, bam, gvcf, vcf -> [ meta, reads ] }, params.skip_trimming )
    ch_versions      = ch_versions.mix(FASTQ_QC_TRIM.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_QC_TRIM.out.multiqc.collect{it[1]})

    //
    // STAGE 1b: Decide GPU vs CPU toolchain (Parabricks when a capable GPU is present)
    //
    CHECK_GPU ()
    ch_gpu_ok   = CHECK_GPU.out.decision.map { it.text.trim() == 'true' }.first()
    ch_versions = ch_versions.mix(CHECK_GPU.out.versions)
    CHECK_GPU.out.report.subscribe { log.info "[genomeportrait] GPU check:\n${it.text.trim()}" }

    //
    // STAGE 2: Align + mark duplicates + small-variant calling (GPU Parabricks or CPU path)
    //
    ALIGN_AND_CALL ( FASTQ_QC_TRIM.out.reads, ch_prealigned, ch_reference, ch_gpu_ok )
    ch_versions      = ch_versions.mix(ALIGN_AND_CALL.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(ALIGN_AND_CALL.out.metrics.collect{it[1]})

    // Statistics for call sets that were supplied rather than made here
    BCFTOOLS_STATS_INPUT ( ch_portrait_vcf )
    ch_versions     = ch_versions.mix(BCFTOOLS_STATS_INPUT.out.versions.first())
    ch_report_parts = ch_report_parts.mix(
        ALIGN_AND_CALL.out.stats.mix(BCFTOOLS_STATS_INPUT.out.stats).map{ m,f -> [m,'variant_stats',f] })

    // Every interpretation stage waits for the report container (built in PREPARE_GENOME),
    // so report tasks never race ahead of the image they run in.
    ch_build_ready = ch_reference.report_sif.ifEmpty('ready').first()

    ch_analysis_bam = gateOn( ALIGN_AND_CALL.out.bam.mix(ch_prealigned, ch_portrait_bam), ch_build_ready )  // [ meta, bam/cram, index ]
    ch_vcf          = gateOn( ALIGN_AND_CALL.out.vcf.mix(ch_portrait_vcf), ch_build_ready )                 // [ meta, vcf.gz, tbi ]
    ch_gvcf         = gateOn( ALIGN_AND_CALL.out.gvcf.mix(ch_portrait_gvcf), ch_build_ready )               // [ meta, g.vcf.gz, tbi ]
    // A real BAM where one exists (Qualimap cannot read CRAM)
    ch_real_bam     = ALIGN_AND_CALL.out.md_bam.mix(
        ch_prealigned.mix(ch_portrait_bam).filter { meta, f, idx -> f.name.endsWith('.bam') } )

    //
    // STAGE 3: Coverage, and what it says about the genome (sex chromosomes, mtDNA copies)
    //
    MOSDEPTH ( ch_analysis_bam, ch_reference.fasta )
    GENOME_FACTS (
        MOSDEPTH.out.summary.join(MOSDEPTH.out.global_dist),
        file("${projectDir}/scripts/genome_facts.py", checkIfExists: true)
    )
    ch_versions      = ch_versions.mix(MOSDEPTH.out.versions.first(), GENOME_FACTS.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.summary.collect{it[1]}, MOSDEPTH.out.global_dist.collect{it[1]})
    ch_report_parts  = ch_report_parts.mix(GENOME_FACTS.out.tsv.map{ m,f -> [m,'genome_facts',f] })
    // [ meta.id, 'XX' | 'XY' | 'unknown' ]: samplesheet value if given, else inferred from depth
    ch_sex = GENOME_FACTS.out.sex.map { meta, f -> [ meta.id, f.text.trim() ] }

    // DNA from a second person in the library: checked before anything is read into the genome
    if (!params.skip_contamination) {
        VERIFYBAMID2 ( ch_analysis_bam, ch_reference.fasta, ch_reference.fai )
        CONTAMINATION_REPORT ( VERIFYBAMID2.out.selfsm, file("${projectDir}/scripts/contamination_report.py", checkIfExists: true) )
        ch_versions     = ch_versions.mix(VERIFYBAMID2.out.versions.first(), CONTAMINATION_REPORT.out.versions.first())
        ch_report_parts = ch_report_parts.mix(CONTAMINATION_REPORT.out.tsv.map{ m,f -> [m,'contamination_facts',f] })
    }

    if (!params.skip_qc) {
        BAM_QC ( ch_analysis_bam, ch_real_bam, ch_reference )
        ch_versions      = ch_versions.mix(BAM_QC.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_QC.out.multiqc.collect{it[1]})
    }

    //
    // STAGE 5: SV / CNV
    //
    if (params.call_sv || params.call_cnv) {
        BAM_CALL_SV_CNV ( ch_analysis_bam, ch_reference )
        ch_versions     = ch_versions.mix(BAM_CALL_SV_CNV.out.versions)
        ch_report_parts = ch_report_parts.mix(BAM_CALL_SV_CNV.out.summary.map{ m,f -> [m,'sv_cnv',f] })
    }

    //
    // STAGE 6: Annotation (VEP + gnomAD/ClinVar/dbNSFP)
    //
    if (!params.skip_annotation) {
        VCF_ANNOTATE ( ch_vcf, ch_reference )
        ch_versions      = ch_versions.mix(VCF_ANNOTATE.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(VCF_ANNOTATE.out.reports.collect{it[1]})
        ch_report_parts  = ch_report_parts.mix(VCF_ANNOTATE.out.tsv.map{ m,f -> [m,'annotation',f] })
    }

    //
    // STAGE 7: Ancestry & admixture
    //
    if (!params.skip_ancestry) {
        BAM_ANCESTRY ( ch_analysis_bam, ch_vcf, ch_reference )
        ch_versions     = ch_versions.mix(BAM_ANCESTRY.out.versions)
        ch_report_parts = ch_report_parts.mix(BAM_ANCESTRY.out.results.map{ m,f -> [m,'ancestry',f] })
    }

    //
    // STAGE 7b: Neanderthal / Denisovan segments (hmmix on the person's own variants)
    //
    if (!params.skip_archaic) {
        ARCHAIC_HMMIX (
            ch_vcf,
            ch_reference.hmmix_resources,
            file("${projectDir}/scripts/run_hmmix.sh", checkIfExists: true),
            file("${projectDir}/scripts/archaic_report.py", checkIfExists: true),
            file("${projectDir}/assets/archaic/1000g_reference.tsv", checkIfExists: true),
            ch_reference.hmmix_ready
        )
        ch_versions     = ch_versions.mix(ARCHAIC_HMMIX.out.versions.first())
        ch_report_parts = ch_report_parts.mix(ARCHAIC_HMMIX.out.results.map{ m,f -> [m,'archaic',f] })
    }

    //
    // STAGE 7c: Chromosome painting (local ancestry: Beagle phasing, RFMix)
    //
    if (!params.skip_painting) {
        PAINT_CHROMOSOMES (
            ch_gvcf,
            ch_reference.painting_reference,
            file("${projectDir}/scripts/pgs_engine.py", checkIfExists: true),
            file("${projectDir}/scripts/run_painting.sh", checkIfExists: true),
            file("${projectDir}/scripts/painting_report.py", checkIfExists: true),
            ch_reference.painting_ready
        )
        ch_versions     = ch_versions.mix(PAINT_CHROMOSOMES.out.versions.first())
        ch_report_parts = ch_report_parts.mix(PAINT_CHROMOSOMES.out.results.map{ m,f -> [m,'painting',f] })
    }

    //
    // STAGE 8: Maternal / paternal lineage haplogroups
    //
    if (!params.skip_haplogroups) {
        BAM_HAPLOGROUPS ( ch_analysis_bam, ch_vcf, ch_reference, ch_sex )
        ch_versions     = ch_versions.mix(BAM_HAPLOGROUPS.out.versions)
        ch_report_parts = ch_report_parts.mix(BAM_HAPLOGROUPS.out.results.map{ m,f -> [m,'haplogroups',f] })
    }

    //
    // STAGE 9: Forensic STR profile (CODIS) + repeat expansions
    //
    if (!params.skip_str) {
        BAM_FORENSIC_STR ( ch_analysis_bam, ch_reference, ch_sex )
        ch_versions     = ch_versions.mix(BAM_FORENSIC_STR.out.versions)
        ch_report_parts = ch_report_parts.mix(BAM_FORENSIC_STR.out.results.map{ m,f -> [m,'forensic_str',f] })
    }

    //
    // STAGE 10: Polygenic scores — ancestry-adjusted, with direction and effect size
    //
    if (!params.skip_prs) {
        VCF_POLYGENIC ( ch_gvcf, ch_reference.panel )
        ch_versions     = ch_versions.mix(VCF_POLYGENIC.out.versions)
        ch_report_parts = ch_report_parts.mix(VCF_POLYGENIC.out.fragment.map{ m,f -> [m,'polygenic',f] })
    }

    //
    // STAGE 10b: the genome as a 23andMe-format raw data file, at the GSA array's SNPs
    //
    if (!params.skip_23andme) {
        ARRAY_PANEL (
            [ params.array_manifest_b37, params.array_manifest_b38 ],
            file("${projectDir}/scripts/array_panel.py", checkIfExists: true)
        )
        EXPORT_23ANDME (
            ch_gvcf.map { meta, g, t -> [ meta.id, meta, g, t ] }
                .join(ch_sex)
                .map { id, meta, g, t, sex -> [ meta, g, t, sex ] },
            ARRAY_PANEL.out.sites,
            ch_reference.fasta, ch_reference.fai,
            file("${projectDir}/scripts/export_23andme.py", checkIfExists: true),
            file("${projectDir}/scripts/pgs_engine.py", checkIfExists: true)
        )
        ch_versions = ch_versions.mix(ARRAY_PANEL.out.versions, EXPORT_23ANDME.out.versions.first())
    }

    //
    // STAGE 10c: blood groups (ABO, Kell, Duffy, Kidd from genotypes; RhD from RHD depth)
    //
    if (!params.skip_blood_groups) {
        RHD_DEPTH ( ch_analysis_bam, ch_reference.fasta, ch_reference.fai,
                    file("${projectDir}/assets/blood_groups/rhd_regions.bed", checkIfExists: true) )
        BLOOD_GROUPS (
            ch_gvcf.map { meta, g, t -> [ meta.id, meta, g, t ] }
                .join( RHD_DEPTH.out.tsv.map { meta, f -> [ meta.id, f ] } )
                .map { id, meta, g, t, rhd -> [ meta, g, t, rhd ] },
            file("${projectDir}/assets/blood_groups/sites.pvar", checkIfExists: true),
            file("${projectDir}/scripts/pgs_engine.py", checkIfExists: true),
            file("${projectDir}/scripts/blood_groups.py", checkIfExists: true)
        )
        ch_versions     = ch_versions.mix(RHD_DEPTH.out.versions.first(), BLOOD_GROUPS.out.versions.first())
        ch_report_parts = ch_report_parts.mix(BLOOD_GROUPS.out.tsv.map{ m,f -> [m,'blood_groups',f] })
    }

    //
    // STAGE 11: HLA typing (T1K)
    //
    if (!params.skip_hla) {
        BAM_HLA ( ch_analysis_bam, ch_reference )
        ch_versions     = ch_versions.mix(BAM_HLA.out.versions)
        ch_report_parts = ch_report_parts.mix(BAM_HLA.out.results.map{ m,f -> [m,'hla',f] })
    }

    //
    // STAGE 11b: KIR gene content (T1K)
    //
    if (!params.skip_kir) {
        BAM_KIR ( ch_analysis_bam, ch_reference )
        ch_versions     = ch_versions.mix(BAM_KIR.out.versions)
        ch_report_parts = ch_report_parts.mix(BAM_KIR.out.results.map{ m,f -> [m,'kir',f] })
    }

    //
    // STAGE 12: Pharmacogenomics (PharmCAT)
    //
    if (!params.skip_pharmcat) {
        BAM_PHARMACOGENOMICS ( ch_analysis_bam, ch_reference )
        ch_versions     = ch_versions.mix(BAM_PHARMACOGENOMICS.out.versions)
        ch_report_parts = ch_report_parts.mix(BAM_PHARMACOGENOMICS.out.results.map{ m,f -> [m,'pharmcat',f] })
    }

    //
    // STAGE 13: "Popular" non-clinical trait lookup
    //
    if (!params.skip_traits) {
        VCF_TRAITS ( ch_analysis_bam, ch_reference, file(params.trait_panel) )
        ch_versions     = ch_versions.mix(VCF_TRAITS.out.versions)
        ch_report_parts = ch_report_parts.mix(VCF_TRAITS.out.results.map{ m,f -> [m,'traits',f] })
    }

    //
    // STAGE 14: Extras (telomere length, mtDNA heteroplasmy)
    //
    BAM_EXTRAS ( ch_analysis_bam, ch_reference )
    ch_versions     = ch_versions.mix(BAM_EXTRAS.out.versions)
    ch_report_parts = ch_report_parts.mix(BAM_EXTRAS.out.results.map{ m,f -> [m,'extras',f] })

    //
    // Collate software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'genomeportrait_software_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // STAGE 15a: MultiQC
    //
    ch_multiqc_report = Channel.empty()
    if (!params.skip_multiqc) {
        ch_multiqc_config = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
        summary_params    = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
        ch_multiqc_files = ch_multiqc_files
            .mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
            .mix(ch_collated_versions)

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            [], [], [], []
        )
        ch_multiqc_report = MULTIQC.out.report
        ch_versions       = ch_versions.mix(MULTIQC.out.versions)
    }

    //
    // STAGE 15b: Unified GenomePortrait HTML report
    //   Group all per-sample report parts into one bundle per sample.
    //
    ch_report_parts
        .map { meta, tag, f -> [ meta, [ tag, f ] ] }
        .groupTuple()
        .map { meta, parts -> [ meta, parts.collect{ it[1] } ] }
        .set { ch_report_bundle }

    GENOMEPORTRAIT_REPORT (
        ch_report_bundle,
        ch_multiqc_report.ifEmpty([]),
        file("${projectDir}/scripts/genomeportrait_report.py", checkIfExists: true)
    )
    ch_versions = ch_versions.mix(GENOMEPORTRAIT_REPORT.out.versions)

    emit:
    multiqc_report = ch_multiqc_report
    portrait       = GENOMEPORTRAIT_REPORT.out.html
    versions       = ch_versions
}
