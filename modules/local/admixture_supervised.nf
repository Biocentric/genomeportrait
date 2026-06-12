process ADMIXTURE_SUPERVISED {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::admixture=1.3.0 bioconda::plink2=2.00a5.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/admixture:1.3.0--0' :
        'biocontainers/admixture:1.3.0--0' }"

    input:
    tuple val(meta), path(pgen), path(pvar), path(psam)
    tuple path(ref_pgen), path(ref_pvar), path(ref_psam)

    output:
    tuple val(meta), path("*.Q.tsv"), emit: q
    path "versions.yml",              emit: versions

    script:
    // Supervised mode: reference individuals carry population labels in a .pop file;
    // the query individual is left unlabelled and its ancestry fractions are estimated.
    """
    plink2 --pfile ${pgen.baseName} --make-bed --out merged
    # Build .pop: superpopulation for reference samples, '-' for the query
    awk 'NR>1 {print (\$1=="${meta.id}" ? "-" : \$5)}' ${ref_psam} > merged.pop
    K=\$(sort -u merged.pop | grep -v '^-\$' | wc -l)
    admixture --supervised -j${task.cpus} merged.bed \$K
    paste <(awk 'NR>1{print \$2}' ${psam}) merged.\${K}.Q > ${meta.id}.Q.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        admixture: \$(admixture --version 2>&1 | grep -oP 'ADMIXTURE Version \\K[0-9.]+')
    END_VERSIONS
    """
}
