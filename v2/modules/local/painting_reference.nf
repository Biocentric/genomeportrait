process PAINTING_REFERENCE {
    tag "1kGP phased panel"
    label 'process_medium'
    label 'process_long'
    storeDir "${params.reference_base}/painting"

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? (params.painting_sif ?: "${params.reference_base}/containers/genomeportrait-painting-${params.painting_image_version}.sif") : params.painting_docker }"

    input:
    path  build_script   // scripts/paint_reference.sh
    path  holdout        // assets/painting/holdout.tsv
    val   ready          // the locally built image exists

    output:
    path "reference",    emit: reference
    path "versions.yml", emit: versions

    script:
    // The 1000 Genomes 30x phased panel (~30 GB of VCF) is streamed one autosome at a time,
    // MD5-checked, and cut to ~900,000 common SNPs: Beagle's phasing reference (all genomes
    // but the hold-outs) and RFMix's (400 unadmixed genomes), ~2 GB in all. See the script.
    """
    bash ${build_script} ${params.painting_panel_url} ${params.painting_ped_url} ${params.painting_map_url} \\
        ${holdout} ${task.cpus} reference

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        panel: 1kGP_high_coverage_Illumina phased SNV_INDEL_SV panel (20220422)
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p "reference"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: true
    END_VERSIONS
    """
}
