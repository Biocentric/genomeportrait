process DOWNLOAD_KGP_HGDP {
    tag "kgp_hgdp_panel"
    label 'process_medium'
    label 'process_long'
    storeDir storedir

    conda "bioconda::plink2=2.00a5.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0' :
        'biocontainers/plink2:2.00a5.10--h4ac6f70_0' }"

    input:
    tuple val(pgen_url), val(pvar_url), val(psam_url)
    val storedir

    output:
    tuple path("hgdp_1kg.pgen"), path("hgdp_1kg.pvar"), path("hgdp_1kg.psam"), emit: panel
    path "versions.yml", emit: versions

    script:
    """
    wget --quiet --continue -O panel.pgen.zst "${pgen_url}"
    wget --quiet --continue -O panel.pvar.zst "${pvar_url}"
    wget --quiet --continue -O hgdp_1kg.psam  "${psam_url}"
    plink2 --zst-decompress panel.pgen.zst > hgdp_1kg.pgen
    plink2 --zst-decompress panel.pvar.zst > hgdp_1kg.pvar

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ .*//')
    END_VERSIONS
    """
}
