process CLAIR3 {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::clair3=1.0.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/clair3:1.0.10--py39hf5e1c6e_0' :
        'biocontainers/clair3:1.0.10--py39hf5e1c6e_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("*.clair3.vcf.gz"), path("*.clair3.vcf.gz.tbi"), emit: vcf
    path "versions.yml",                                                   emit: versions

    script:
    // Long-read model selection; assumes Clair3 models bundled in container
    def platform = meta.platform == 'pacbio' ? 'hifi' : 'ont'
    """
    MODEL=\$(find / -type d -name "${platform}*" -path "*models*" 2>/dev/null | head -n1)
    run_clair3.sh \\
        --bam_fn=$bam --ref_fn=$fasta \\
        --threads=$task.cpus --platform=$platform \\
        --model_path=\${MODEL} --output=clair3_out
    cp clair3_out/merge_output.vcf.gz ${meta.id}.clair3.vcf.gz
    cp clair3_out/merge_output.vcf.gz.tbi ${meta.id}.clair3.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(run_clair3.sh --version 2>&1 | sed 's/Clair3 v//')
    END_VERSIONS
    """
}
