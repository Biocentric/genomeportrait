process ADMIXTURE_SUPERVISED {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::admixture=1.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/admixture:1.3.0--0' :
        'biocontainers/admixture:1.3.0--0' }"

    input:
    tuple val(meta), path(c_bed), path(c_bim), path(c_fam)
    tuple path(panel_bed), path(panel_bim), path(panel_fam)
    path  ref_psam

    output:
    tuple val(meta), path("${meta.id}.admixture.tsv"), emit: q
    path "versions.yml",                               emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Supervised: reference rows carry their population label (.pop, in combined.fam order),
    // the query is '-'. ADMIXTURE doesn't name its Q columns, so emit per-individual
    // IID + label + Q1..QK and let the report map column->population and pick the query
    // (which the report finds by IID, since plink/admixture row order != input order).
    // psam cols: #IID PAT MAT SEX SuperPop(5) Population(6).
    def lcol = params.ancestry_label_col == 'Population' ? 6 : 5
    """
    # ADMIXTURE reads <bedprefix>.pop, so the .pop MUST be named to match the bed file.
    awk 'NR==FNR{lab[\$1]=\$${lcol}; next}{print ((\$2 in lab)? lab[\$2] : "-")}' \\
        <(tail -n +2 ${ref_psam}) ${meta.id}.combined.fam > ${meta.id}.combined.pop
    K=\$(grep -vx - ${meta.id}.combined.pop | sort -u | wc -l)
    echo "ADMIXTURE K=\$K (${params.ancestry_label_col})"

    admixture --supervised -j${task.cpus} ${meta.id}.combined.bed \$K

    # Assemble IID + label + Q with a single awk (the admixture container lacks paste/cat/tr);
    # all three files share combined.fam row order, so join by line number.
    awk -v K=\$K '
      BEGIN { printf "IID\\tlabel"; for (j=1;j<=K;j++) printf "\\tQ%d", j; printf "\\n" }
      FNR==1 { fi++ }
      fi==1 { iid[FNR]=\$2; next }
      fi==2 { pop[FNR]=\$1; next }
      fi==3 { printf "%s\\t%s", iid[FNR], pop[FNR]; for (j=1;j<=NF;j++) printf "\\t%s", \$j; printf "\\n" }
    ' ${meta.id}.combined.fam ${meta.id}.combined.pop ${meta.id}.combined.\${K}.Q > ${meta.id}.admixture.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        admixture: \$(admixture --version 2>&1 | grep -oP 'ADMIXTURE Version \\K[0-9.]+')
    END_VERSIONS
    """
}
