process BCFTOOLS_GENOTYPE_SITES {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bcftools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.21--h8b25389_0' :
        'biocontainers/bcftools:1.21--h8b25389_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fai
    path  sites          // CHROM POS REF ALT (chr-prefixed)

    output:
    tuple val(meta), path("${meta.id}.anc.vcf.gz"), path("${meta.id}.anc.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Genotype the sample at the ancestry SNPs straight from the alignment, for the panel's
    // own alleles (-C alleles; -A makes it write the given ALT at hom-ref sites). Without that, a hom-ref site is called with ALT ".", its ID
    // (chr:pos:ref:.) never matches the panel, and it is dropped at the merge: in v1 30,312
    // of 50,000 sites were lost that way, leaving mostly variant sites and skewing admixture
    // and kinship. The sample is renamed to the samplesheet name: the read group's SM may be
    // a different (real) name.
    """
    cut -f1,2 $sites > regions.tsv
    sort -k1,1 -k2,2n $sites | bgzip -c > alleles.tsv.gz
    tabix -s1 -b2 -e2 alleles.tsv.gz
    bcftools mpileup -R regions.tsv -f $fasta -a FORMAT/AD,FORMAT/DP -Ou $bam 2>/dev/null \\
        | bcftools call -m -A -C alleles -T alleles.tsv.gz -Oz -o raw.vcf.gz
    echo "${meta.id}" > name.txt
    bcftools reheader -s name.txt -o ${meta.id}.anc.vcf.gz raw.vcf.gz
    bcftools index -t ${meta.id}.anc.vcf.gz
    rm -f raw.vcf.gz
    echo "genotyped \$(bcftools view -H ${meta.id}.anc.vcf.gz | wc -l) sites; without a panel allele: \$(bcftools view -H ${meta.id}.anc.vcf.gz | awk '\$5 == "."' | wc -l)" >&2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    """
    echo "" | gzip > "${meta.id}.anc.vcf.gz"
    touch "${meta.id}.anc.vcf.gz.tbi"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: true
    END_VERSIONS
    """
}
