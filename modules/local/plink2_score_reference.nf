process PLINK2_SCORE_REFERENCE {
    tag "pgs_reference"
    label 'process_medium'
    label 'process_long'
    // The reference distribution depends only on the PGS files, never on the sample, so
    // compute it once and reuse it for every genome. Keyed by the requested PGS IDs.
    storeDir "${params.reference_base}/pgs_reference/${params.pgs_ids.replaceAll(/[^A-Za-z0-9]/, '_')}"

    conda "bioconda::plink2=2.00a5.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0' :
        'biocontainers/plink2:2.00a5.10--h4ac6f70_0' }"

    input:
    tuple path(ref_pgen), path(ref_pvar), path(ref_psam)
    path  scorefile

    output:
    path "ref_scores",   emit: scores
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // A raw polygenic score is just a weighted allele sum: its magnitude tracks how many
    // variants the score contains, so on its own it means nothing. Scoring the 1000G+HGDP
    // reference panel with the SAME files gives a distribution to place the sample in.
    // Subset the panel to the PGS variants FIRST (one pass over the large panel), then score
    // each PGS against that small fileset.
    """
    mkdir -p ref_scores
    zcat -f $scorefile > combined.txt

    # union of scored positions, as <chr>:<pos> to match --set-all-var-ids '@:#'
    awk -F'\\t' 'NR>1 {print \$1":"\$2}' combined.txt | sort -T . -u > pgs_union.ids
    echo "PGS variants (union): \$(wc -l < pgs_union.ids)" >&2

    plink2 --pfile ${ref_pgen.baseName} --allow-extra-chr \\
        --set-all-var-ids '@:#' --rm-dup force-first \\
        --extract pgs_union.ids --make-pgen --out panel_pgs
    echo "panel subset: \$(grep -vc '^#' panel_pgs.pvar) variants" >&2

    for acc in \$(tail -n +2 combined.txt | cut -f8 | sort -u); do
        safe=\$(echo "\$acc" | tr -c 'A-Za-z0-9._-' '_')
        { printf "ID\\tA1\\tBETA\\n"; awk -F'\\t' -v a="\$acc" 'NR>1 && \$8==a {print \$1":"\$2"\\t"\$3"\\t"\$5}' combined.txt; } > sf.\$safe.txt
        plink2 --pfile panel_pgs \\
            --score sf.\$safe.txt 1 2 3 header-read no-mean-imputation cols=+scoresums,+denom ignore-dup-ids \\
            --out ref_scores/\$safe 2>/dev/null || echo "reference scoring failed for \$acc" >&2
    done
    cp ${ref_psam} ref_scores/reference.psam
    ls ref_scores | awk 'NR<=20' >&2 || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ .*//')
    END_VERSIONS
    """
}
