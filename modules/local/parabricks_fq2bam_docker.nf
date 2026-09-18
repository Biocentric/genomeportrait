process PARABRICKS_FQ2BAM_DOCKER {
    tag "$meta.id"
    label 'process_high'
    label 'process_gpu'
    // NOT stageInMode 'copy': the trimmed FASTQ pair for a WGS sample is ~95 GB and
    // copying it into the task dir would double that on disk for no benefit. Inputs
    // stay symlinked, and params.parabricks_docker_binds mounts the underlying
    // filesystems at their real paths so the symlink targets resolve inside Docker.

    // Deliberately NOT a Nextflow-managed container process.
    //
    // Parabricks ignores CUDA_VISIBLE_DEVICES: it always takes CUDA device 0.
    // Singularity's `--nv` exposes every GPU on the host, so under Singularity
    // pbrun lands on whichever card is first on the PCI bus — which on a mixed
    // box may be a small display GPU, and fq2bam then dies with
    //   cudaSafeCall() failed at samGenerator.cu: out of memory
    // Singularity can only isolate a single GPU through `--nvccli`, which needs
    // `nvidia-container-cli` configured in the root-owned singularity.conf.
    //
    // `docker run --gpus device=N` does that isolation with no root setup, so
    // this module runs natively and shells out to Docker. Everything else in the
    // pipeline keeps using the configured container engine. Requires the invoking
    // user to be in the `docker` group.

    input:
    tuple val(meta), path(reads), path(interval_file)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(index)            // classic BWA index dir (.amb/.ann/.bwt/.pac/.sa)
    tuple val(meta4), path(known_sites)      // optional VCF(s) for the BQSR table

    output:
    tuple val(meta), path("*.bam")                  , emit: bam,               optional:true
    tuple val(meta), path("*.bai")                  , emit: bai,               optional:true
    tuple val(meta), path("*.table")                , emit: bqsr_table,        optional:true
    tuple val(meta), path("*.duplicate-metrics.txt"), emit: duplicate_metrics, optional:true
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                  = task.ext.args ?: ''
    def prefix                = task.ext.prefix ?: "${meta.id}"
    def in_fq_command         = meta.single_end ? "--in-se-fq ${reads}" : "--in-fq ${reads}"
    def known_sites_command   = known_sites ? known_sites.collect { "--knownSites $it" }.join(' ') + " --out-recal-file ${prefix}.table" : ''
    def interval_file_command = interval_file ? interval_file.collect { "--interval-file $it" }.join(' ') : ''
    // pbrun defaults the read group to LB:lib1 PL:bar. PL:bar is not a valid platform and
    // GATK rejects it; set the same @RG the CPU path writes via bwa-mem2 -R.
    def platform              = (meta.platform ?: 'illumina').toString().toUpperCase()
    def num_gpus              = task.accelerator ? "--num-gpus ${task.accelerator.request}" : '--num-gpus 1'
    def image                 = params.parabricks_container
    def device                = params.parabricks_gpu_device
    // Mount the underlying filesystems at their real paths so symlinked inputs
    // (reads, reference, index) resolve to the same absolute path inside Docker.
    def binds                 = (params.parabricks_docker_binds ?: '').toString()
    // Docker defaults starve pbrun when it is not run as root: /dev/shm is 64 MB and
    // locked memory is capped at 8 MB. Root can raise memlock (CAP_IPC_LOCK); an
    // unprivileged uid cannot, so CUDA pinned-host allocations fail as the dataset
    // grows -- silently, producing a header-only BAM.
    def dopts                 = (params.parabricks_docker_opts ?: '').toString()
    """
    # pbrun takes the BWA index prefix from --ref, i.e. it opens
    # <ref>.{amb,ann,bwt,pac,sa} beside the reference. BWA_INDEX emits them into a
    # bwa/ subdirectory, so move them next to the staged FASTA first.
    # stageInMode is 'copy', so these are task-local files.
    INDEX=\$(find -L . -name "*.amb" | head -n1 | sed 's/\\.amb\$//')
    if [ -n "\${INDEX}" ] && [ "\${INDEX#./}" != "${fasta}" ]; then
        for e in amb ann bwt pac sa; do
            # Symlink, never move: BWA_INDEX's output is staged as a symlink to its
            # storeDir, so `mv` here follows it and rips the index out of the
            # reference store, silently arming a 90-minute rebuild on the next run.
            [ -e "${fasta}.\${e}" ] || ln -sf "\$(readlink -f "\${INDEX}.\${e}")" "${fasta}.\${e}"
        done
    fi

    mkdir -p pbtmp

    docker run --rm --gpus device=${device} \\
        -u \$(id -u):\$(id -g) ${dopts} \\
        ${binds} \\
        -v "\$PWD":"\$PWD" -w "\$PWD" \\
        ${image} \\
        pbrun fq2bam \\
            --ref ${fasta} \\
            ${in_fq_command} \\
            --read-group-sm ${meta.id} \\
            --read-group-lb ${meta.id} \\
            --read-group-pl ${platform} \\
            --read-group-id-prefix ${meta.id} \\
            ${known_sites_command} \\
            ${interval_file_command} \\
            --out-bam "\$PWD/${prefix}.bam" \\
            --tmp-dir "\$PWD/pbtmp" \\
            ${num_gpus} \\
            $args

    # pbrun can exit 0 having written nothing but a BAM header. Observed on a
    # 1.24 B-read WGS run: every phase reported success, <prefix>_chrs.txt tallied
    # 1.1 B mapped reads per chromosome, and the BAM held 0 records. Nextflow saw a
    # zero exit and passed the empty file downstream, where it cost hours before
    # anything noticed. Check the output against pbrun's own read tally and fail
    # loudly instead, keeping --tmp-dir for diagnosis.
    bam_bytes=\$(stat -c %s "${prefix}.bam")
    if [ -s "${prefix}_chrs.txt" ]; then
        mapped=\$(awk '{s+=\$2} END{print s+0}' "${prefix}_chrs.txt")
        min_bytes=\$(( mapped * 10 ))        # real BAMs run ~80-120 bytes/read
        if [ "\$bam_bytes" -lt "\$min_bytes" ]; then
            echo "ERROR: pbrun exited 0 but ${prefix}.bam is \$bam_bytes bytes for \$mapped mapped reads." >&2
            echo "       Expected at least \$min_bytes. The BAM is header-only - refusing to continue." >&2
            echo "       Temp dir kept at \$PWD/pbtmp for inspection." >&2
            exit 65
        fi
    fi

    rm -rf pbtmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbrun: \$(docker run --rm ${image} pbrun version 2>&1 | grep -oP 'pbrun: \\K[0-9.\\-]+' | head -n1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bam ${prefix}.bam.bai
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pbrun: stub
    END_VERSIONS
    """
}
