process PHARMCAT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::pharmcat=2.15.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://pgkb/pharmcat:2.15.5' :
        'pgkb/pharmcat:2.15.5' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("*.report.json"), emit: json
    tuple val(meta), path("*.report.html"), emit: html
    path "versions.yml",                    emit: versions

    script:
    """
    # Genotype the sample AT PharmCAT's own positions straight from the alignment, WITHOUT -v
    # so homozygous-reference sites are emitted. Feeding the DeepVariant VCF instead is
    # useless here: it is variant-sites-only, so every PGx position where the sample is
    # hom-ref is simply absent and PharmCAT cannot tell "reference" from "no data" — it found
    # 13 of 1132 positions and called nearly every gene Unknown/No Result.
    POS=/pharmcat/pharmcat_positions.vcf
    bcftools query -f '%CHROM\\t%POS\\n' \$POS | sort -k1,1 -k2,2n -u > pgx_sites.tsv
    echo "genotyping \$(wc -l < pgx_sites.tsv) PharmCAT positions from the alignment" >&2

    bcftools mpileup -R pgx_sites.tsv -f $fasta -a FORMAT/DP -Ou $bam 2>/dev/null \\
        | bcftools call -m -Oz -o pgx_raw.vcf.gz

    # Strip AD/PL. AD is Number=R and hom-ref records carry ALT="." (so R=1), which makes
    # PharmCAT's normaliser abort with "could not merge FORMAT tag AD" and preprocess 0
    # records. PharmCAT only needs GT; DP is harmless and kept.
    bcftools annotate -x FORMAT/AD,FORMAT/PL -Oz -o ${meta.id}.pgx.vcf.gz pgx_raw.vcf.gz
    bcftools index -t ${meta.id}.pgx.vcf.gz

    # Preprocess to PharmCAT-ready VCF, then run the full pipeline.
    # (pharmcat_pipeline uses its bundled reference positions — no --reference flag.)
    # -reporterJson is REQUIRED: without it PharmCAT writes only the HTML report, the
    # .report.json copy below falls back to "{}", and the downstream reporter then dies with
    # "No columns to parse from file".
    pharmcat_pipeline ${meta.id}.pgx.vcf.gz -o pharmcat_out -bf ${meta.id} -reporterJson || true

    # Fall back to phenotype.json (always written) if the reporter JSON is missing.
    if [ -s pharmcat_out/${meta.id}.report.json ]; then
        cp pharmcat_out/${meta.id}.report.json ${meta.id}.report.json
    elif [ -s pharmcat_out/${meta.id}.phenotype.json ]; then
        echo "no reporter JSON; falling back to phenotype.json" >&2
        cp pharmcat_out/${meta.id}.phenotype.json ${meta.id}.report.json
    else
        cp pharmcat_out/${meta.id}*.report.json ${meta.id}.report.json 2>/dev/null || echo '{}' > ${meta.id}.report.json
    fi
    cp pharmcat_out/${meta.id}*.report.html ${meta.id}.report.html 2>/dev/null || echo '<html><body>PharmCAT produced no report</body></html>' > ${meta.id}.report.html
    rm -f pgx_raw.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pharmcat: \$(echo 2.15.5)
        bcftools: \$(bcftools --version | head -n1 | sed 's/bcftools //')
    END_VERSIONS
    """
}
