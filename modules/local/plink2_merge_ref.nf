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
    // The panel ships rsID variant IDs while the sample VCF has none, so re-ID BOTH to
    // <chr>:<pos> (plink2 normalises chr1->1, so the schemes line up across GRCh38) and
    // intersect on that. --allow-extra-chr tolerates the alt/random contigs in the sample
    // VCF; --autosome keeps 1-22 for PCA/ADMIXTURE.
    """
    plink2 --vcf $vcf --allow-extra-chr --autosome \\
        --max-alleles 2 --min-alleles 2 --snps-only \\
        --set-all-var-ids '@:#' --rm-dup force-first \\
        --make-pgen --out sample

    plink2 --pfile ${ref_pgen.baseName} --allow-extra-chr --autosome \\
        --max-alleles 2 --min-alleles 2 --snps-only \\
        --set-all-var-ids '@:#' --rm-dup force-first \\
        --make-pgen --out ref

    # variant IDs (chr:pos) shared by sample and panel
    awk '!/^#/{print \$3}' sample.pvar | sort -T . -u > sample.ids
    awk '!/^#/{print \$3}' ref.pvar    | sort -T . -u > ref.ids
    comm -12 sample.ids ref.ids > common.ids
    echo "shared variants: \$(wc -l < common.ids)"

    plink2 --pfile ref    --extract common.ids --make-pgen --out ref_common
    plink2 --pfile sample --extract common.ids --make-pgen --out sample_common
    plink2 --pfile ref_common --pmerge sample_common --make-pgen --out ${meta.id}.merged

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ .*//')
    END_VERSIONS
    """
}
