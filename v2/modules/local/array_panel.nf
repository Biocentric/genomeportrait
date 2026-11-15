process ARRAY_PANEL {
    tag "GSA-24 v3.0"
    label 'process_single'
    // Sample-independent: the export's site table is built once from Illumina's manifests
    // (downloaded here, not shipped with the pipeline) and cached.
    storeDir "${params.reference_base}/array_panel/gsa_v3_r2"   // r2: deterministic choice among duplicate probes

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? (params.report_sif ?: "${params.reference_base}/containers/genomeportrait-report-${params.report_image_version}.sif") : params.report_docker }"

    input:
    val   urls      // [ GRCh37 manifest zip, GRCh38 manifest zip ]
    path  script

    output:
    path "gsa_v3_sites.tsv.gz", emit: sites
    path "versions.yml",        emit: versions

    script:
    """
    python3 - <<'PY'
import urllib.request, zipfile
for name, url in (("b37.zip", "${urls[0]}"), ("b38.zip", "${urls[1]}")):
    urllib.request.urlretrieve(url, name)
    with zipfile.ZipFile(name) as z:
        z.extractall(name[:-4])
PY
    python3 ${script} --b37 \$(ls b37/*.csv | head -1) --b38 \$(ls b38/*.csv | head -1) --out gsa_v3_sites.tsv.gz
    rm -rf b37 b38 b37.zip b38.zip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        manifest: GSA-24v3-0_A1 (GRCh37), GSA-24v3-0_A2 (GRCh38)
    END_VERSIONS
    """

    stub:
    """
    echo "" | gzip > "gsa_v3_sites.tsv.gz"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: true
    END_VERSIONS
    """
}
