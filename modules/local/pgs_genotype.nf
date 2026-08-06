process PGS_GENOTYPE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bcftools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.21--h8b25389_0' :
        'biocontainers/bcftools:1.21--h8b25389_0' }"

    input:
    tuple val(meta), path(gvcf), path(tbi)
    path  fasta
    path  fai
    path  scorefile

    output:
    tuple val(meta), path("${meta.id}.pgs.vcf.gz"), path("${meta.id}.pgs.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // A PGS is a weighted sum over FIXED positions, so every scored position must have a
    // genotype. The DeepVariant VCF holds variant sites only: wherever the sample is
    // homozygous-reference the record is simply absent, and plink2 reads that as missing
    // rather than as a dosage. Because PGS effect alleles are frequently the REFERENCE
    // allele (where hom-ref means a maximal dosage of 2), the loss is directional and drags
    // the score down — on this sample only ~50% of each score's positions were present and
    // every score fell below the entire reference distribution.
    //
    // The gVCF already encodes hom-ref as blocks, so expand it instead of re-piling the BAM:
    // restricting to the scored blocks FIRST and expanding only those took ~5 s per 31k
    // positions, versus ~6 min for the equivalent bcftools mpileup (~3 min vs ~4 h genome-wide).
    """
    CHRPFX=\$(bcftools view -h ${gvcf} | grep -m1 '^##contig=<ID=chr' >/dev/null 2>&1 && echo chr || echo '')
    zcat -f ${scorefile} | awk -F'\\t' -v p="\$CHRPFX" 'NR>1 {print p\$1"\\t"\$2}' \\
        | sort -k1,1 -k2,2n -u > pgs_sites.tsv
    echo "PGS positions requested: \$(wc -l < pgs_sites.tsv)" >&2

    bcftools view -R pgs_sites.tsv -Ou ${gvcf} 2>/dev/null \\
        | bcftools convert --gvcf2vcf -f ${fasta} 2>/dev/null \\
        | bcftools view -T pgs_sites.tsv -Oz -o ${meta.id}.pgs.vcf.gz 2>/dev/null
    bcftools index -t ${meta.id}.pgs.vcf.gz
    echo "genotyped: \$(bcftools view -H ${meta.id}.pgs.vcf.gz | wc -l)" >&2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
