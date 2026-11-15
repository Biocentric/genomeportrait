process BUILD_PAINTING_CONTAINER {
    tag "genomeportrait-painting:${params.painting_image_version}"
    label 'process_low'
    label 'process_long'

    // Runs on the HOST (no container) — invokes apptainer/singularity/docker. Beagle and RFMix
    // are on bioconda, but not in one image with curl and Python; this builds one. storeDir caches it.
    storeDir "${params.reference_base}/containers"

    input:
    path def_file        // containers/painting/Singularity.def

    output:
    path "genomeportrait-painting-${params.painting_image_version}.sif", emit: sif, optional: true
    path "painting_container_${params.painting_image_version}.ready",    emit: ready
    path "versions.yml",                                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ver = params.painting_image_version
    def sif = "genomeportrait-painting-${ver}.sif"
    def engine = workflow.containerEngine ?: 'singularity'
    if (engine in ['singularity', 'apptainer'])
        """
        BUILDER=\$(command -v apptainer || command -v singularity)
        if [ -z "\$BUILDER" ]; then
            echo "ERROR: neither apptainer nor singularity found on PATH" >&2
            exit 1
        fi
        echo "Building ${sif} with \$BUILDER (Beagle, RFMix, bcftools) ..." >&2
        \$BUILDER build --fakeroot ${sif} ${def_file} \\
            || \$BUILDER build ${sif} ${def_file}
        touch painting_container_${ver}.ready

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            builder: \$(\$BUILDER --version 2>&1 | head -n1)
        END_VERSIONS
        """
    else if (engine == 'docker')
        """
        echo "Building docker image ${params.painting_docker} ..." >&2
        docker build -t ${params.painting_docker} \$(dirname \$(readlink -f ${def_file}))
        touch painting_container_${ver}.ready

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            builder: \$(docker --version)
        END_VERSIONS
        """
    else
        """
        echo "Container engine '${engine}': the painting image is not packaged as one biocontainer — build with singularity/docker." >&2
        touch painting_container_${ver}.ready
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            builder: "none (${engine})"
        END_VERSIONS
        """

    stub:
    // the .ready marker only: an empty image in reference_base/containers would be taken for a real one
    """
    touch "painting_container_${params.painting_image_version}.ready"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: true
    END_VERSIONS
    """
}
