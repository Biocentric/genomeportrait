process BUILD_HIPSTR_CONTAINER {
    tag "genomeportrait-hipstr:${params.hipstr_image_version}"
    label 'process_low'
    label 'process_long'

    // Runs on the HOST (no container) — invokes apptainer/singularity/docker to build
    // HipSTR from source (not on bioconda). storeDir caches it across runs.
    storeDir "${params.reference_base}/containers"

    input:
    path def_file        // containers/hipstr/Singularity.def

    output:
    path "genomeportrait-hipstr-${params.hipstr_image_version}.sif", emit: sif, optional: true
    path "hipstr_container_${params.hipstr_image_version}.ready",    emit: ready
    path "versions.yml",                                             emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ver = params.hipstr_image_version
    def sif = "genomeportrait-hipstr-${ver}.sif"
    def engine = workflow.containerEngine ?: 'singularity'
    if (engine in ['singularity', 'apptainer'])
        """
        BUILDER=\$(command -v apptainer || command -v singularity)
        if [ -z "\$BUILDER" ]; then
            echo "ERROR: neither apptainer nor singularity found on PATH" >&2
            exit 1
        fi
        echo "Building ${sif} with \$BUILDER (compiles HipSTR from source) ..." >&2
        \$BUILDER build --fakeroot ${sif} ${def_file} \\
            || \$BUILDER build ${sif} ${def_file}
        touch hipstr_container_${ver}.ready

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            builder: \$(\$BUILDER --version 2>&1 | head -n1)
        END_VERSIONS
        """
    else if (engine == 'docker')
        """
        echo "Building docker image ${params.hipstr_docker} ..." >&2
        docker build -t ${params.hipstr_docker} \$(dirname \$(readlink -f ${def_file}))
        touch hipstr_container_${ver}.ready

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            builder: \$(docker --version)
        END_VERSIONS
        """
    else
        """
        echo "Container engine '${engine}': HipSTR has no conda package — build with singularity/docker." >&2
        touch hipstr_container_${ver}.ready
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            builder: "none (${engine})"
        END_VERSIONS
        """
}
