process CYRIUS {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? (params.cyrius_sif ?: "${params.reference_base}/containers/genomeportrait-cyrius-${params.cyrius_image_version}.sif") : params.cyrius_docker }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fai
    val   ready      // the locally built image exists

    output:
    tuple val(meta), path("${meta.id}.cyrius.tsv"),       emit: tsv
    tuple val(meta), path("${meta.id}.cyrius.json"),      emit: json, optional: true
    tuple val(meta), path("${meta.id}.cyp2d6_outside.tsv"), emit: outside
    path "versions.yml",                                  emit: versions

    script:
    // Cyrius reads the CYP2D6/CYP2D7 locus plus ~3,000 normalisation regions for depth.
    // A PASS call becomes a PharmCAT outside call (gene <tab> diplotype). Ambiguous calls
    // ("a;b") and no-calls are not passed on: PharmCAT then reports CYP2D6 as not called
    // instead of guessing.
    """
    echo "\$PWD/${bam}" > manifest.txt
    cyrius --manifest manifest.txt --genome 38 --prefix ${meta.id} --outDir cyrius \\
        --threads ${task.cpus} --reference ${fasta}
    cp cyrius/${meta.id}.tsv ${meta.id}.cyrius.tsv
    [ -f cyrius/${meta.id}.json ] && cp cyrius/${meta.id}.json ${meta.id}.cyrius.json || true
    awk -F'\\t' 'NR == 2 && \$3 == "PASS" && \$2 != "None" && \$2 !~ /;/ { print "CYP2D6\\t" \$2 }' \\
        ${meta.id}.cyrius.tsv > ${meta.id}.cyp2d6_outside.tsv
    echo "Cyrius: \$(awk -F'\\t' 'NR == 2 { print \$2, \$3 }' ${meta.id}.cyrius.tsv)" >&2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cyrius: 1.1.1
    END_VERSIONS
    """
}
