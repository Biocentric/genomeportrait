process PLINK_ANCESTRY_MERGE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::plink=1.90b6.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b7.7--h18e278d_1' :
        'biocontainers/plink:1.90b7.7--h18e278d_1' }"

    input:
    tuple val(meta), path(s_bed), path(s_bim), path(s_fam)
    tuple path(panel_bed), path(panel_bim), path(panel_fam)

    output:
    tuple val(meta), path("${meta.id}.combined.bed"), path("${meta.id}.combined.bim"), path("${meta.id}.combined.fam"), emit: combined
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // plink2 --pmerge is unimplemented, so merge with plink 1.9 --bmerge. Force the sample's
    // A1 to the panel coding first; on the 3+-allele error, exclude the missnp and re-merge.
    def panel = panel_bed.baseName
    """
    plink --bfile ${meta.id}.sample --a1-allele ${panel_bim} 5 2 --make-bed --allow-no-sex --out s_al

    if ! plink --bfile ${panel} --bmerge s_al --make-bed --allow-no-sex --out ${meta.id}.combined 2>/dev/null; then
        plink --bfile ${panel} --exclude ${meta.id}.combined-merge.missnp --make-bed --allow-no-sex --out ref_x
        plink --bfile s_al      --exclude ${meta.id}.combined-merge.missnp --make-bed --allow-no-sex --out s_x
        plink --bfile ref_x --bmerge s_x --make-bed --allow-no-sex --out ${meta.id}.combined
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(plink --version 2>&1 | head -n1 | grep -oE 'v[0-9.]+[a-z0-9.-]*' | head -1)
    END_VERSIONS
    """
}
