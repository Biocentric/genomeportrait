process PLINK2_ROH {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::plink=1.90b7.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b7.2--h4ac6f70_0' :
        'biocontainers/plink:1.90b7.2--h4ac6f70_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.roh.tsv"), emit: roh
    path "versions.yml",                emit: versions

    script:
    def args = task.ext.args ?: '--homozyg'
    """
    plink --vcf $vcf --double-id $args --out ${meta.id}
    # Inbreeding coefficient F (method-of-moments) + ROH burden
    plink --vcf $vcf --double-id --het --out ${meta.id}.het
    {
        echo -e "metric\\tvalue"
        echo -e "n_roh_segments\\t\$(tail -n +2 ${meta.id}.hom 2>/dev/null | wc -l)"
        echo -e "total_roh_kb\\t\$(awk 'NR>1{s+=\$9}END{print s+0}' ${meta.id}.hom 2>/dev/null)"
        echo -e "inbreeding_F\\t\$(awk 'NR>1{print \$6}' ${meta.id}.het.het 2>/dev/null | head -n1)"
    } > ${meta.id}.roh.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(plink --version 2>&1 | sed 's/^PLINK v//; s/ .*//')
    END_VERSIONS
    """
}
