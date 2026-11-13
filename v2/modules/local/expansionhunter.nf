process EXPANSIONHUNTER {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::expansionhunter=5.0.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/expansionhunter:5.0.0--hf366f20_0' :
        'biocontainers/expansionhunter:5.0.0--hf366f20_0' }"

    input:
    tuple val(meta), path(bam), path(bai), val(sex_chromosomes)
    path  fasta
    path  fai
    path  catalog

    output:
    tuple val(meta), path("*.expansionhunter.vcf"), emit: vcf
    path "versions.yml",                            emit: versions

    script:
    // ExpansionHunter's --sex only sets the ploidy it assumes on chrX/chrY loci. v1 took it
    // from the samplesheet and fell back to 'female' for an undeclared sex, which genotypes
    // an XY sample's X-linked repeats as diploid. v2 passes the sex resolved from depth.
    def sex = sex_chromosomes == 'XY' ? 'male' : 'female'
    """
    ExpansionHunter \\
        --reads $bam \\
        --reference $fasta \\
        --variant-catalog $catalog \\
        --sex $sex \\
        --output-prefix ${meta.id}.expansionhunter

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        expansionhunter: \$(ExpansionHunter --version 2>&1 | grep -oP '[0-9.]+' | head -n1)
    END_VERSIONS
    """
}
