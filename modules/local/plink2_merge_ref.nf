process PLINK2_MERGE_REF {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::plink2=2.00a5.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0' :
        'biocontainers/plink2:2.00a5.10--h4ac6f70_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple path(ref_pgen), path(ref_pvar), path(ref_psam)

    output:
    tuple val(meta), path("${meta.id}.merged.pgen"), path("${meta.id}.merged.pvar"), path("${meta.id}.merged.psam"), emit: merged
    path "versions.yml", emit: versions

    script:
    """
    # Convert sample VCF -> pgen, restrict to bi-allelic SNPs shared with the panel
    plink2 --vcf $vcf --max-alleles 2 --snps-only --set-all-var-ids '@:#' \\
        --make-pgen --out sample

    # Intersect on variant IDs and merge sample into the reference panel
    plink2 --pfile ${ref_pgen.baseName} --extract sample.pvar --make-pgen --out ref_common
    plink2 --pfile sample --extract ref_common.pvar --make-pgen --out sample_common
    plink2 --pfile ref_common --pmerge sample_common --make-pgen --out ${meta.id}.merged

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ .*//')
    END_VERSIONS
    """
}
