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
    awk 'NR==FNR{lab[\$1]=\$${lcol}; next}{print ((\$2 in lab)? lab[\$2] : "-")}' \\
        <(tail -n +2 ${ref_psam}) ${meta.id}.combined.fam > combined.pop
    K=\$(grep -vx - combined.pop | sort -u | wc -l)
    echo "ADMIXTURE K=\$K (${params.ancestry_label_col})"

    admixture --supervised -j${task.cpus} ${meta.id}.combined.bed \$K

    { printf "IID\\tlabel"; for i in \$(seq 1 \$K); do printf "\\tQ%d" \$i; done; printf "\\n";
      paste <(awk '{print \$2}' ${meta.id}.combined.fam) combined.pop ${meta.id}.combined.\${K}.Q | tr ' ' '\\t';
    } > ${meta.id}.admixture.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        admixture: \$(admixture --version 2>&1 | grep -oP 'ADMIXTURE Version \\K[0-9.]+')
    END_VERSIONS
    """
}
