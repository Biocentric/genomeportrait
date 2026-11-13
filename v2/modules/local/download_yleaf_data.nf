process DOWNLOAD_YLEAF_DATA {
    tag "yleaf_data"
    label 'process_single'
    storeDir "${params.reference_base}/yleaf_data"

    conda "bioconda::yleaf=3.2.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/yleaf:3.2.1--pyh1286868_0' :
        'biocontainers/yleaf:3.2.1--pyh1286868_0' }"

    output:
    path "yleaf_data", emit: data
    path "versions.yml", emit: versions

    script:
    // The yleaf biocontainer ships only the .py files — its `data/` directory (Y-SNP position
    // tables, ISOGG/FTDNA/YFull marker sets and the haplogroup prediction tree) is missing, so
    // Yleaf crashes before it can call anything. Fetch that directory from the upstream repo
    // once and cache it; it is small (~200 KB).
    """
    # Plain extract then move: the container's BusyBox tar has no --wildcards.
    wget -q -O yleaf_src.tar.gz "${params.yleaf_data_url}"
    mkdir -p src && tar -xzf yleaf_src.tar.gz -C src
    mv "\$(find src -type d -path '*/yleaf/data' | head -1)" yleaf_data
    rm -rf src yleaf_src.tar.gz
    echo "yleaf data files: \$(find yleaf_data -type f | wc -l)"
    test -f yleaf_data/hg38/new_positions.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yleaf_data: ${params.yleaf_data_url}
    END_VERSIONS
    """
}
