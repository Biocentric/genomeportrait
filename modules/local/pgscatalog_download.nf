process PGSCATALOG_DOWNLOAD {
    tag "$pgs_ids"
    label 'process_single'

    conda "bioconda::pgscatalog-utils=1.4.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pgscatalog-utils:1.4.4--pyhdfd78af_0' :
        'biocontainers/pgscatalog-utils:1.4.4--pyhdfd78af_0' }"

    input:
    val pgs_ids

    output:
    path "scorefiles_combined.txt.gz", emit: scorefiles
    path "pgs_metadata.tsv",           emit: metadata
    path "versions.yml",               emit: versions

    script:
    """
    # pgscatalog-utils 1.4.x requires the output dirs to already exist
    mkdir -p scorefiles meta

    # Download scoring files (already on GRCh38, so no liftover needed)
    pgscatalog-download --pgs ${pgs_ids.replace(',', ' ')} \\
        --build GRCh38 --outdir scorefiles

    # Combine into a single weighted-allele file for PLINK2 --score
    pgscatalog-combine -s scorefiles/*.txt.gz \\
        --target-build GRCh38 \\
        -o scorefiles_combined.txt.gz

    # Lightweight metadata table for the report (trait labels per PGS ID)
    pgscatalog-download --pgs ${pgs_ids.replace(',', ' ')} --metadata-only -o meta || true
    printf "pgs_id\\ttrait\\n" > pgs_metadata.tsv
    for id in \$(echo ${pgs_ids} | tr ',' ' '); do
        printf "%s\\t%s\\n" "\$id" "\$(grep -h "\$id" scorefiles/*.txt 2>/dev/null | grep -m1 trait_reported | cut -f2- || echo 'see PGS Catalog')" >> pgs_metadata.tsv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pgscatalog-utils: \$(pgscatalog-download --version 2>&1 | grep -oP '[0-9.]+' | head -n1)
    END_VERSIONS
    """
}
