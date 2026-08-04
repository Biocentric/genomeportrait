process BUILD_YLEAF_CONTAINER {
    tag "genomeportrait-yleaf:${params.yleaf_image_version}"
    label 'process_low'
    label 'process_long'

    // Runs on the HOST (no container) — invokes apptainer/singularity/docker. The published
    // yleaf biocontainer has no data/ dir and no samtools/bcftools, so it cannot run at all;
    // this builds a usable one. storeDir caches it across runs.
    storeDir "${params.reference_base}/containers"

    input:
    path def_file        // containers/yleaf/Singularity.def

    output:
    path "genomeportrait-yleaf-${params.yleaf_image_version}.sif", emit: sif, optional: true
    path "yleaf_container_${params.yleaf_image_version}.ready",    emit: ready
    path "versions.yml",                                           emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ver = params.yleaf_image_version
    def sif = "genomeportrait-yleaf-${ver}.sif"
    def engine = workflow.containerEngine ?: 'singularity'
    if (engine in ['singularity', 'apptainer'])
        """
        BUILDER=\$(command -v apptainer || command -v singularity)
        if [ -z "\$BUILDER" ]; then
            echo "ERROR: neither apptainer nor singularity found on PATH" >&2
            exit 1
        fi
        echo "Building ${sif} with \$BUILDER (yleaf + samtools + bcftools) ..." >&2
        \$BUILDER build --fakeroot ${sif} ${def_file} \\
            || \$BUILDER build ${sif} ${def_file}
        touch yleaf_container_${ver}.ready

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            builder: \$(\$BUILDER --version 2>&1 | head -n1)
        END_VERSIONS
        """
    else if (engine == 'docker')
        """
        echo "Building docker image ${params.yleaf_docker} ..." >&2
        docker build -t ${params.yleaf_docker} \$(dirname \$(readlink -f ${def_file}))
        touch yleaf_container_${ver}.ready

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            builder: \$(docker --version)
        END_VERSIONS
        """
    else
        """
        echo "Container engine '${engine}': the yleaf biocontainer is unusable — build with singularity/docker." >&2
        touch yleaf_container_${ver}.ready
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            builder: "none (${engine})"
        END_VERSIONS
        """
}
