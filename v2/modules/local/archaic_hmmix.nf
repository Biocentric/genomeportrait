process ARCHAIC_HMMIX {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? (params.hmmix_sif ?: "${params.reference_base}/containers/genomeportrait-hmmix-${params.hmmix_image_version}.sif") : params.hmmix_docker }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple path(mask), path(mutrates), path(outgroup), path(ancestral), path(archaic)
    path  run_script     // scripts/run_hmmix.sh
    path  report_script  // scripts/archaic_report.py
    path  reference      // assets/archaic/1000g_reference.tsv
    val   ready          // the locally built image exists

    output:
    tuple val(meta), path("${meta.id}.archaic"),                  emit: results
    tuple val(meta), path("${meta.id}.archaic_segments.raw.tsv"), emit: segments
    tuple val(meta), path("${meta.id}.trained.json"),             emit: params
    path "versions.yml",                                          emit: versions

    script:
    // Archaic (Neanderthal / Denisovan) segments: scripts/run_hmmix.sh runs hmmix
    // create_ingroup (per autosome, in parallel), train and decode on the person's own
    // variants; archaic_report.py turns the segments into the report section.
    """
    mkdir -p res
    for f in ${mask} ${mutrates} ${outgroup} ${ancestral} ${archaic}; do ln -s "\$(readlink -f \$f)" res/; done
    bash ${run_script} ${vcf} res ${task.cpus} ${meta.id}
    python3 ${report_script} --sample ${meta.id} --decoded ${meta.id}.decoded.diploid.txt --mask ${mask} \\
        --outdir ${meta.id}.archaic --segments-out ${meta.id}.archaic_segments.raw.tsv --reference ${reference}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmix: \$(pip show hmmix 2>/dev/null | sed -n 's/^Version: //p')
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p "${meta.id}.archaic"
    touch "${meta.id}.archaic_segments.raw.tsv"
    touch "${meta.id}.trained.json"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: true
    END_VERSIONS
    """
}
