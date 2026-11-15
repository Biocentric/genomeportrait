process PAINT_CHROMOSOMES {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? (params.painting_sif ?: "${params.reference_base}/containers/genomeportrait-painting-${params.painting_image_version}.sif") : params.painting_docker }"

    input:
    tuple val(meta), path(gvcf), path(tbi)
    path  reference      // PAINTING_REFERENCE
    path  engine         // scripts/pgs_engine.py: genotypes at fixed sites from the gVCF, hom-ref included
    path  run_script     // scripts/run_painting.sh
    path  report_script  // scripts/painting_report.py
    val   ready          // the locally built image exists

    output:
    tuple val(meta), path("${meta.id}.painting"),                  emit: results
    tuple val(meta), path("${meta.id}.painting_segments.raw.tsv"), emit: segments
    tuple val(meta), path("rfmix/*"),                              emit: rfmix
    path "versions.yml",                                           emit: versions

    script:
    // The person's genotypes at the reference sites come from the gVCF, so hom-ref sites
    // count (a variant-only VCF would leave most sites missing); Beagle phases them against
    // the 1000 Genomes panel, RFMix paints each phased copy, painting_report.py summarises.
    """
    python3 ${engine} target --gvcf ${gvcf} --pvar ${reference}/sites.pvar --sample ${meta.id} \\
        --threads ${task.cpus} --out target.vcf.gz --stats target_stats.tsv
    mkdir -p rfmix
    bash ${run_script} target.vcf.gz ${reference} ${task.cpus} rfmix/${meta.id}
    python3 ${report_script} --sample ${meta.id} --msp-prefix rfmix/${meta.id} --outdir ${meta.id}.painting \\
        --segments-out ${meta.id}.painting_segments.raw.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        beagle: 5.4_22Jul22.46e
        rfmix: 2.03-r0
    END_VERSIONS
    """
}
