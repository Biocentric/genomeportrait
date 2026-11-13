process BUILD_REPORT_CONTAINER {
    tag "genomeportrait-report:${params.report_image_version}"
    label 'process_low'

    // Runs on the HOST (no container) because it invokes apptainer/singularity/docker.
    // storeDir caches the built image so it is only built once, then re-used by every run.
    storeDir "${params.reference_base}/containers"

    input:
    path def_file        // containers/report/Singularity.def
    path env_file        // containers/report/environment.yml

    output:
    path "genomeportrait-report-${params.report_image_version}.sif",  emit: sif, optional: true
    path "report_container_${params.report_image_version}.ready",     emit: ready
    path "versions.yml",                                              emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ver = params.report_image_version
    def sif = "genomeportrait-report-${ver}.sif"
    def engine = workflow.containerEngine ?: 'singularity'
    if (engine in ['singularity', 'apptainer'])
        """
        # Prefer apptainer (unprivileged --fakeroot builds), fall back to singularity.
        BUILDER=\$(command -v apptainer || command -v singularity)
        if [ -z "\$BUILDER" ]; then
            echo "ERROR: neither apptainer nor singularity found on PATH" >&2
            exit 1
        fi
        echo "Building ${sif} with \$BUILDER ..." >&2
        \$BUILDER build --fakeroot ${sif} ${def_file} \\
            || \$BUILDER build ${sif} ${def_file}
        touch report_container_${ver}.ready

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            builder: \$(\$BUILDER --version 2>&1 | head -n1)
        END_VERSIONS
        """
    else if (engine == 'docker')
        """
        echo "Building docker image ${params.report_docker} ..." >&2
        docker build -t ${params.report_docker} \$(dirname \$(readlink -f ${def_file}))
        touch report_container_${ver}.ready

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            builder: \$(docker --version)
        END_VERSIONS
        """
    else
        // conda / mamba profile: nothing to build, the conda directive provides the env
        """
        echo "Container engine '${engine}': using conda env, no image build needed." >&2
        touch report_container_${ver}.ready
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            builder: "none (conda profile)"
        END_VERSIONS
        """
}
