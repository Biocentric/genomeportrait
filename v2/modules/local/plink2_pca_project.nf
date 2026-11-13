process PLINK2_PCA_PROJECT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::plink2=2.00a5.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0' :
        'biocontainers/plink2:2.00a5.10--h4ac6f70_0' }"

    input:
    tuple val(meta), path(c_bed), path(c_bim), path(c_fam)
    tuple path(panel_bed), path(panel_bim), path(panel_fam)

    output:
    tuple val(meta), path("${meta.id}.pca.tsv"), emit: eigenvec
    path "versions.yml",                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Learn the PCA basis on the reference panel, then PROJECT everyone (reference + query)
    // onto it via --score (allele-orientation safe). The .eigenvec.allele file has a
    // PROVISIONAL_REF? column, so A1=col6 and the PC weights are cols 7-16.
    def panel = panel_bed.baseName
    """
    plink2 --bfile ${panel} --freq --out reff
    plink2 --bfile ${panel} --read-freq reff.afreq --pca 10 allele-wts --out refpca

    plink2 --bfile ${meta.id}.combined --read-freq reff.afreq \\
        --score refpca.eigenvec.allele 2 6 header-read no-mean-imputation variance-standardize \\
        --score-col-nums 7-16 --out proj

    # tidy: IID PC1..PC10 (drop FID/ALLELE_CT/dosage cols; PCs are the last 10 columns)
    awk 'NR==1{printf "IID"; for(i=1;i<=10;i++) printf "\\tPC%d", i; printf "\\n"; next}
         {printf "%s", \$2; for(i=NF-9;i<=NF;i++) printf "\\t%s", \$i; printf "\\n"}' proj.sscore > ${meta.id}.pca.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ .*//')
    END_VERSIONS
    """
}
