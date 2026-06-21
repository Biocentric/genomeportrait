process HIPSTR {
    tag "$meta.id"
    label 'process_medium'

    // HipSTR is not on bioconda — use the image built at runtime by BUILD_HIPSTR_CONTAINER.
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? params.hipstr_sif : params.hipstr_docker }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fai
    path  str_bed
    val   ready          // gate: ensures BUILD_HIPSTR_CONTAINER finished before this runs

    output:
    tuple val(meta), path("*.hipstr.vcf.gz"), emit: vcf
    path "versions.yml",                      emit: versions

    script:
    def panel = params.str_panel ?: 'curated'
    def codis = 'CSF1PO|FGA|TH01|TPOX|vWA|D3S1358|D5S818|D7S820|D8S1179|D13S317|D16S539|D18S51|D21S11|D1S1656|D2S441|D2S1338|D10S1248|D12S391|D19S433|D22S1045'
    """
    # Build the STR regions for HipSTR. The reference bed is genome-wide (~1.6M loci, hours);
    #   curated = named forensic/genealogy markers (CODIS + Y-STRs + Marshfield, ~860) [default]
    #   codis   = CODIS core only (~20)
    #   all     = full genome-wide reference
    BED=${str_bed}
    if [[ \$BED == *.gz ]]; then gunzip -kf \$BED; BED=${str_bed.baseName}; fi
    case "${panel}" in
        all)   REGIONS=\$BED ;;
        codis) grep -iwE '${codis}' \$BED > str_regions.bed; REGIONS=str_regions.bed ;;
        *)     grep -vE 'Human_STR_[0-9]' \$BED > str_regions.bed; REGIONS=str_regions.bed ;;
    esac
    echo "HipSTR panel='${panel}': \$(wc -l < \$REGIONS) loci" >&2

    HipSTR \\
        --bams $bam \\
        --fasta $fasta \\
        --regions \${REGIONS} \\
        --min-reads ${params.hipstr_min_reads} \\
        --str-vcf ${meta.id}.hipstr.vcf.gz \\
        --output-filters

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hipstr: \$(HipSTR --version 2>&1 | grep -oP '[0-9.]+' | head -n1)
    END_VERSIONS
    """
}
