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
        -t GRCh38 \\
        -o scorefiles_combined.txt.gz

    # Lightweight metadata table for the report (trait labels per PGS ID).
    # NB the downloaded scoring files are .txt.gz — the previous glob looked for *.txt, found
    # nothing, and every trait fell through to the "see PGS Catalog" placeholder.
    pgscatalog-download --pgs ${pgs_ids.replace(',', ' ')} --metadata-only -o meta || true
    printf "pgs_id\ttrait\tn_variants\turl\n" > pgs_metadata.tsv
    for id in \$(echo ${pgs_ids} | tr ',' ' '); do
        f=\$(ls scorefiles/\${id}*.txt.gz 2>/dev/null | head -1)
        trait="."; nvar="."
        if [ -n "\$f" ]; then
            hdr=\$(zcat "\$f" | head -40)
            trait=\$(printf '%s' "\$hdr" | grep -m1 '^#trait_reported=' | cut -d= -f2- || true)
            [ -z "\$trait" ] && trait=\$(printf '%s' "\$hdr" | grep -m1 '^#trait_mapped=' | cut -d= -f2- || true)
            nvar=\$(printf '%s' "\$hdr" | grep -m1 '^#variants_number=' | cut -d= -f2- || true)
        fi
        [ -z "\$trait" ] && trait="(trait not stated in scoring file)"
        [ -z "\$nvar" ] && nvar="."
        printf "%s\t%s\t%s\thttps://www.pgscatalog.org/score/%s/\n" "\$id" "\$trait" "\$nvar" "\$id" >> pgs_metadata.tsv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pgscatalog-utils: \$(pgscatalog-download --version 2>&1 | grep -oP '[0-9.]+' | head -n1)
    END_VERSIONS
    """
}
