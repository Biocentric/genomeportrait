process PGS_PANEL_DOWNLOAD {
    tag "HGDP+1kGP"
    label 'process_low'
    label 'process_long'
    storeDir "${params.reference_base}/pgs/reference_panel"

    conda "bioconda::plink2=2.00a5.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0' :
        'biocontainers/plink2:2.00a5.10--h4ac6f70_0' }"

    input:
    val url

    output:
    tuple path("GRCh38_HGDP+1kGP_ALL.pgen"), path("GRCh38_HGDP+1kGP_ALL.pvar.zst"),
          path("GRCh38_HGDP+1kGP_ALL.psam"), path("GRCh38_HGDP+1kGP.king.cutoff.out.id"), emit: panel
    path "versions.yml",                                                                  emit: versions

    script:
    // The archive (16 GB) holds the panel for GRCh37 AND GRCh38, ~30 GB unpacked. Stream it
    // through plink2's zstd decoder into tar and keep only the GRCh38 members, so neither
    // an intermediate .tar nor the unused build ever lands on disk.
    // A truncated download fails deep inside a run, not at the start: check the size.
    """
    wget --quiet --continue -O panel.tar.zst "${url}"
    size=\$(wc -c < panel.tar.zst)
    if [ "\$size" -lt 15000000000 ]; then
        echo "reference panel download is \$size bytes; expected ~15.95 GB" >&2
        exit 1
    fi
    plink2 --zst-decompress panel.tar.zst | tar -xf - \\
        'GRCh38_HGDP+1kGP_ALL.pgen' 'GRCh38_HGDP+1kGP_ALL.pvar.zst' \\
        'GRCh38_HGDP+1kGP_ALL.psam' 'GRCh38_HGDP+1kGP.king.cutoff.out.id'
    rm -f panel.tar.zst

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ .*//')
        panel: pgsc_HGDP+1kGP_v1
    END_VERSIONS
    """
}
