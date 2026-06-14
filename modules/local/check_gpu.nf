process CHECK_GPU {
    tag "${params.use_gpu}"
    label 'process_single'
    // Runs on the host (no container) so it can query the host nvidia-smi.

    output:
    path "gpu.decision",  emit: decision   // contains "true" or "false"
    path "gpu.report.txt", emit: report
    path "versions.yml",  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def mode    = params.use_gpu          // 'off' | 'on' | 'auto'
    def min_mem = params.gpu_min_memory_mb
    def min_cc  = params.gpu_min_compute_cap
    """
    set +e
    MODE="${mode}"
    REPORT=gpu.report.txt
    : > \$REPORT

    if [ "\$MODE" = "off" ]; then
        echo "use_gpu=off -> CPU tools" | tee -a \$REPORT
        echo false > gpu.decision
    else
        if ! command -v nvidia-smi >/dev/null 2>&1; then
            echo "nvidia-smi not found." | tee -a \$REPORT
            if [ "\$MODE" = "on" ]; then echo "ERROR: --use_gpu on but no GPU runtime." | tee -a \$REPORT; exit 64; fi
            echo false > gpu.decision
        else
            nvidia-smi --query-gpu=name,memory.total,compute_cap --format=csv,noheader >> \$REPORT 2>/dev/null
            # Capable if ANY GPU has >= min_mem MiB AND compute capability >= min_cc
            CAPABLE=\$(nvidia-smi --query-gpu=memory.total,compute_cap --format=csv,noheader,nounits 2>/dev/null | \\
                awk -F', ' -v m=${min_mem} -v c=${min_cc} 'BEGIN{ok=0} {gsub(/ /,"",\$1); if (\$1+0 >= m && \$2+0 >= c) ok=1} END{print ok}')
            if [ "\$CAPABLE" = "1" ]; then
                echo "Parabricks-capable GPU found (>=${min_mem} MiB, compute>=${min_cc}) -> GPU tools" | tee -a \$REPORT
                echo true > gpu.decision
            else
                echo "GPU present but below Parabricks spec (need >=${min_mem} MiB & compute>=${min_cc})." | tee -a \$REPORT
                if [ "\$MODE" = "on" ]; then echo "ERROR: --use_gpu on but no capable GPU." | tee -a \$REPORT; exit 64; fi
                echo "Falling back to CPU tools." | tee -a \$REPORT
                echo false > gpu.decision
            fi
        fi
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nvidia-smi: \$(command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi --version 2>/dev/null | grep -oP 'DRIVER version: \\K[0-9.]+' | head -n1 || echo "n/a")
    END_VERSIONS
    """
}
