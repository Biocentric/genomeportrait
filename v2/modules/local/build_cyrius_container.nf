process BUILD_CYRIUS_CONTAINER {
    tag "genomeportrait-cyrius:${params.cyrius_image_version}"
    label 'process_low'
    label 'process_long'

    // Runs on the HOST (no container) — invokes apptainer/singularity/docker. Cyrius is on
    // PyPI only, so there is no published image; this builds one. storeDir caches it.
    storeDir "${params.reference_base}/containers"

    input:
    path def_file        // containers/cyrius/Singularity.def

    output:
    path "genomeportrait-cyrius-${params.cyrius_image_version}.sif", emit: sif, optional: true
    path "cyrius_container_${params.cyrius_image_version}.ready",    emit: ready
    path "versions.yml",                                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ver = params.cyrius_image_version
    def sif = "genomeportrait-cyrius-${ver}.sif"
    def engine = workflow.containerEngine ?: 'singularity'
    if (engine in ['singularity', 'apptainer'])
        """
        BUILDER=\$(command -v apptainer || command -v singularity)
        if [ -z "\$BUILDER" ]; then
            echo "ERROR: neither apptainer nor singularity found on PATH" >&2
            exit 1
        fi
        echo "Building ${sif} with \$BUILDER (Cyrius + pysam) ..." >&2
        \$BUILDER build --fakeroot ${sif} ${def_file} \\
            || \$BUILDER build ${sif} ${def_file}
        touch cyrius_container_${ver}.ready

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            builder: \$(\$BUILDER --version 2>&1 | head -n1)
        END_VERSIONS
        """
    else if (engine == 'docker')
        """
        echo "Building docker image ${params.cyrius_docker} ..." >&2
        docker build -t ${params.cyrius_docker} \$(dirname \$(readlink -f ${def_file}))
        touch cyrius_container_${ver}.ready

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            builder: \$(docker --version)
        END_VERSIONS
        """
    else
        """
        echo "Container engine '${engine}': Cyrius is not packaged as a biocontainer — build with singularity/docker." >&2
        touch cyrius_container_${ver}.ready
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            builder: "none (${engine})"
        END_VERSIONS
        """
}
