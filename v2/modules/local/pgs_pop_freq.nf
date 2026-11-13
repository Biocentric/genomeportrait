process PGS_POP_FREQ {
    tag "$panel_key"
    label 'process_medium'
    // Allele frequencies of each reference super-population. A sample's missing genotypes are
    // filled from its most similar population's frequencies rather than the pooled panel's.
    storeDir "${params.reference_base}/pgs/${panel_key}/pop_freq"

    conda "bioconda::plink2=2.00a5.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0' :
        'biocontainers/plink2:2.00a5.10--h4ac6f70_0' }"

    input:
    val   panel_key
    tuple path(ref_pgen), path(ref_pvar), path(ref_psam)

    output:
    path "pop_afreq",    emit: afreq
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p pop_afreq
    col=\$(head -1 ${ref_psam} | tr '\\t' '\\n' | grep -nx '${params.pgs_population_col}' | cut -d: -f1)
    for pop in \$(tail -n +2 ${ref_psam} | cut -f\$col | sort -u); do
        awk -F'\\t' -v c=\$col -v p=\$pop 'NR > 1 && \$c == p { print \$1 }' ${ref_psam} > \$pop.ids
        plink2 --threads ${task.cpus} --pfile ${ref_pgen.baseName} --keep \$pop.ids --freq --out pop_afreq/\$pop
    done
    rm -f pop_afreq/*.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ .*//')
    END_VERSIONS
    """
}
