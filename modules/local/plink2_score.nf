process PLINK2_SCORE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::plink2=2.00a5.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0' :
        'biocontainers/plink2:2.00a5.10--h4ac6f70_0' }"

    input:
    tuple val(meta), path(pgen), path(pvar), path(psam)
    path  scorefile

    output:
    tuple val(meta), path("*.sscore"), emit: sscore
    path "versions.yml",               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // pgscatalog-combine long format: chr_name pos effect_allele other_allele effect_weight
    //                                 effect_type is_duplicated accession row_nr
    // Score each PGS separately. The pgen var IDs are <chr>:<pos> (PLINK2_VCF2PGEN, @:#) using
    // the genome's own contig naming, and PGS chr_name uses the same convention, so build the
    // scorefile ID the same way (no chr-prefix munging). no-mean-imputation: with a single
    // sample plink2 cannot impute allele frequencies, so missing genotypes contribute 0.
    """
    zcat -f $scorefile > combined.txt
    for acc in \$(tail -n +2 combined.txt | cut -f8 | sort -u); do
        safe=\$(echo "\$acc" | tr -c 'A-Za-z0-9._-' '_')
        { printf "ID\\tA1\\tBETA\\n"; awk -F'\\t' -v a="\$acc" 'NR>1 && \$8==a {print \$1":"\$2"\\t"\$3"\\t"\$5}' combined.txt; } > sf.\$safe.txt
        plink2 --pfile ${pgen.baseName} \\
            --score sf.\$safe.txt 1 2 3 header-read no-mean-imputation cols=+scoresums,+denom ignore-dup-ids list-variants \\
            --out ${meta.id}.\$safe 2>/dev/null || echo "score failed for \$acc" >&2
    done
    ls *.sscore >/dev/null 2>&1 || printf "" > ${meta.id}.none.sscore

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ .*//')
    END_VERSIONS
    """
}
