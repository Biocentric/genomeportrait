process KING_KINSHIP {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::king=2.3.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/king:2.3.2--h3be2455_0' :
        'biocontainers/king:2.3.2--h3be2455_0' }"

    input:
    tuple val(meta), path(c_bed), path(c_bim), path(c_fam)

    output:
    tuple val(meta), path("*.kin.tsv"), emit: kin
    path "versions.yml",                emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // plink writes FID=0 for everyone, so KING saw all 3200+ reference individuals as ONE
    // family: it put every pair in <prefix>.kin (the WITHIN-family estimator, which assumes
    // the pairs really are relatives and is therefore inflated by shared ancestry) and never
    // produced .kin0. Give each individual its own FID so KING uses the robust BETWEEN-family
    // estimator, which is the one designed for structured/unrelated samples.
    """
    awk '{ \$1 = \$2; print }' OFS='\\t' ${c_fam} > king.fam
    cp ${c_bed} king.bed
    cp ${c_bim} king.bim

    king -b king.bed --kinship --prefix ${meta.id} || true

    # KING emits .kin0 (across families) and/or .kin (within family) — take whichever exists.
    SRC=""
    [ -s ${meta.id}.kin0 ] && SRC=${meta.id}.kin0
    [ -z "\$SRC" ] && [ -s ${meta.id}.kin ] && SRC=${meta.id}.kin

    if [ -n "\$SRC" ]; then
        # Keep pairs involving the query, strongest first, and label the relationship using
        # KING's standard degree cut-offs.
        awk -v s="${meta.id}" -F'\\t' '
          NR==1 { for (i=1;i<=NF;i++) col[\$i]=i; next }
          (\$col["ID1"]==s || \$col["ID2"]==s) {
              k = \$col["Kinship"] + 0
              other = (\$col["ID1"]==s) ? \$col["ID2"] : \$col["ID1"]
              rel = (k>0.354) ? "duplicate/MZ" : (k>0.177) ? "1st-degree" : \\
                    (k>0.0884) ? "2nd-degree" : (k>0.0442) ? "3rd-degree" : "unrelated"
              printf "%s\\t%s\\t%.5f\\t%s\\t%s\\n", s, other, k, rel, \$col["N_SNP"]
          }' "\$SRC" | sort -k3,3gr > pairs.tsv

        {
          printf "sample\\trelative\\tkinship\\tinferred_relationship\\tn_snp\\n"
          if [ -s pairs.tsv ]; then
              head -n 20 pairs.tsv
          else
              printf "%s\\tNA\\tNA\\tno_pairs_computed\\tNA\\n" "${meta.id}"
          fi
        } > ${meta.id}.kin.tsv

        TOP=\$(head -n1 pairs.tsv 2>/dev/null | cut -f3)
        echo "top kinship vs the reference panel: \${TOP:-none} (estimator: \$SRC)" >&2
    else
        printf "sample\\trelative\\tkinship\\tinferred_relationship\\tn_snp\\n%s\\tNA\\tNA\\tking_produced_no_output\\tNA\\n" "${meta.id}" > ${meta.id}.kin.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        king: \$(king 2>&1 | grep -oE 'KING [0-9.]+' | head -n1 | sed 's/KING //')
    END_VERSIONS
    """
}
