process HMMIX_RESOURCES {
    tag "hmmix hg38"
    label 'process_low'
    label 'process_long'
    storeDir "${params.reference_base}/hmmix"

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? (params.hmmix_sif ?: "${params.reference_base}/containers/genomeportrait-hmmix-${params.hmmix_image_version}.sif") : params.hmmix_docker }"

    input:
    path  fetch_script   // scripts/zenodo_fetch.py
    val   ready          // the locally built image exists

    output:
    tuple path("hg38_strick_callability_mask.bed"), path("hg38_mutationrate.bed"),
          path("outgroup"), path("ancestral"), path("archaic"),                       emit: resources
    path "versions.yml",                                                           emit: versions

    script:
    // Skov's hmmix helper files (Zenodo ${params.hmmix_helper_record}: callability mask, mutation
    // rate per Mb, variants seen in the 1000G/HGDP outgroup, Ensembl ancestral alleles) and the
    // four high-coverage archaic genomes (Zenodo ${params.hmmix_archaic_record}: Altai, Vindija
    // and Chagyrskaya Neanderthals, Denisova), ~3.3 GB in all, MD5-checked. The archaic BCFs
    // name chromosomes 1..22; hmmix looks them up as chr1..chr22, so they are renamed here.
    // The outgroup table (1.6 GB) is cut to chrom/pos per chromosome: hmmix hands it to
    // vcftools --exclude-positions, which reads all of it for every call, and the per-sample
    // step runs one call per chromosome.
    """
    python3 ${fetch_script} --record ${params.hmmix_helper_record} --out . \\
        hg38_strick_callability_mask.bed hg38_mutationrate.bed hg38_Outgroup_1000g_HGDP.txt hg38_ancestral.tar.gz
    mkdir -p ancestral
    tar xzf hg38_ancestral.tar.gz -C ancestral --strip-components 1
    rm hg38_ancestral.tar.gz
    test -s ancestral/homo_sapiens_ancestor_22.fa
    mkdir -p outgroup
    awk -F'\t' 'NR > 1 { print \$1 "\t" \$2 > ("outgroup/" \$1 ".txt") }' hg38_Outgroup_1000g_HGDP.txt
    rm hg38_Outgroup_1000g_HGDP.txt
    test -s outgroup/chr22.txt

    python3 ${fetch_script} --record ${params.hmmix_archaic_record} --out archaic_raw \\
        --pattern 'individuals_highcov\\.([0-9]|1[0-9]|2[0-2])\\.bcf'
    for c in \$(seq 1 22); do echo "\$c chr\$c"; done > chr_names.txt
    mkdir -p archaic
    for c in \$(seq 1 22); do
        bcftools annotate --rename-chrs chr_names.txt -Ob -o archaic/individuals_highcov.chr\$c.bcf \\
            --write-index=csi archaic_raw/individuals_highcov.\$c.bcf
    done
    rm -rf archaic_raw chr_names.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmix_helper_files: zenodo.${params.hmmix_helper_record}
        archaic_genomes: zenodo.${params.hmmix_archaic_record}
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //')
    END_VERSIONS
    """

    stub:
    """
    touch "hg38_strick_callability_mask.bed"
    touch "hg38_mutationrate.bed"
    mkdir -p "outgroup"
    mkdir -p "ancestral"
    mkdir -p "archaic"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: true
    END_VERSIONS
    """
}
