process VCF_TO_TABLE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bcftools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.21--h8b25389_0' :
        'biocontainers/bcftools:1.21--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${meta.id}.annotation.tsv.gz"),          emit: tsv
    tuple val(meta), path("${meta.id}.annotation_full.raw.tsv.gz"), emit: full
    path "versions.yml",                                            emit: versions

    script:
    // CLNSIG is added to INFO by bcftools annotate (ClinVar) after VEP, so it is not a CSQ
    // subfield; reference it explicitly as %INFO/CLNSIG.
    """
    # Call-quality columns come last so the existing column positions stay put.
    # DeepVariant emits no MQ/MAPQ field at all, so mapping quality is simply not available
    # here; FILTER + QUAL + GQ + DP + VAF are what it does provide, and together they cover
    # the same "should I believe this call" judgement.
    bcftools +split-vep -d -f \\
      '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%SYMBOL\\t%Consequence\\t%IMPACT\\t%gnomADg_AF\\t%INFO/CLNSIG\\t%Existing_variation\\t[%GT]\\t%FILTER\\t%QUAL\\t[%GQ]\\t[%DP]\\t[%VAF]\\n' \\
      -A tab $vcf 2>/dev/null \\
      | sed '1i CHROM\\tPOS\\tID\\tREF\\tALT\\tGENE\\tCONSEQUENCE\\tIMPACT\\tGNOMAD_AF\\tCLINSIG\\tDBSNP\\tGT\\tFILTER\\tQUAL\\tGQ\\tDEPTH\\tVAF' \\
      | bgzip -c > ${meta.id}.annotation_full.raw.tsv.gz

    # A whole genome yields ~10 million consequence rows, ~99% of them MODIFIER (intergenic /
    # up- and downstream noise). Loading that into the report costs GBs of RAM and buries
    # anything interpretable, so keep the full table as a published artifact and pass on only
    # the rows that carry signal: coding/regulatory impact, or a ClinVar annotation.
    # Two further requirements, both essential:
    #   FILTER=PASS  — ~38% of DeepVariant records are RefCall, the caller concluding the site
    #                  is actually reference.
    #   GT carries an ALT (\$12 ~ /1/) — otherwise the table lists positions where the sample
    #                  is homozygous REFERENCE (0/0) or uncalled (./.). Those are annotated
    #                  sites, not variants this person has: without this, most of the top
    #                  ClinVar "Pathogenic" hits were 0/0, i.e. variants they do NOT carry.
    zcat ${meta.id}.annotation_full.raw.tsv.gz \\
      | awk -F'\\t' 'NR==1 || (\$13=="PASS" && \$12 ~ /1/ && (\$8=="HIGH" || \$8=="MODERATE" || \$8=="LOW" || (\$10!="." && \$10!="")))' \\
      | bgzip -c > ${meta.id}.annotation.tsv.gz
    echo "notable rows: \$(zcat ${meta.id}.annotation.tsv.gz | wc -l) of \$(zcat ${meta.id}.annotation_full.raw.tsv.gz | wc -l)" >&2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
