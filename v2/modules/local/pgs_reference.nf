process PGS_REFERENCE {
    tag "$panel_key"
    label 'process_high'
    storeDir "${params.reference_base}/pgs/${panel_key}/reference"

    conda "bioconda::plink2=2.00a5.10"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.10--h4ac6f70_0' :
        'biocontainers/plink2:2.00a5.10--h4ac6f70_0' }"

    input:
    val   panel_key
    tuple path(pgen), path(pvar_zst), path(psam), path(related)
    path  sites
    path  pca_exclude

    output:
    tuple path("ref.pgen"), path("ref.pvar"), path("ref.psam"),                emit: pfile
    path "ref.afreq",                                                          emit: afreq
    tuple path("pca.acount"), path("pca.eigenvec.allele"), path("pca.cols"),   emit: pca
    path "ref_pcs.sscore",                                                     emit: ref_pcs
    path "versions.yml",                                                       emit: versions

    script:
    def pfx = pgen.baseName
    def mem = task.memory ? "--memory ${(task.memory.toMega() * 0.8) as long}" : ''
    """
    # 1. Reference panel restricted to scored positions, relatives removed (KING cutoff list
    #    shipped with the panel), autosomes only.
    plink2 --threads ${task.cpus} ${mem} --pfile ${pfx} vzs \\
        --extract bed1 ${sites} --remove ${related} --autosome \\
        --make-pgen --out ref
    plink2 --threads ${task.cpus} ${mem} --pfile ref --freq --out ref

    # 2. PCA on common, LD-pruned SNPs outside the long-range LD regions.
    plink2 --threads ${task.cpus} ${mem} --pfile ref \\
        --maf 0.05 --geno 0.01 --snps-only just-acgt --max-alleles 2 \\
        --exclude bed1 ${pca_exclude} --indep-pairwise 1000kb 0.1 --out prune
    plink2 --threads ${task.cpus} ${mem} --pfile ref --extract prune.prune.in \\
        --freq counts --pca ${params.pgs_n_pcs} allele-wts --out pca

    # 3. Project the reference with the SAME command later used for each sample. plink2's
    #    .eigenvec is on a different scale from a --score projection, so mixing the two
    #    would put the sample in the wrong place.
    A1=\$(awk 'NR==1{for(i=1;i<=NF;i++) if(\$i=="A1") print i}' pca.eigenvec.allele)
    P1=\$(awk 'NR==1{for(i=1;i<=NF;i++) if(\$i=="PC1") print i}' pca.eigenvec.allele)
    P2=\$((P1 + ${params.pgs_n_pcs} - 1))
    echo "\$A1 \$P1 \$P2" > pca.cols
    plink2 --threads ${task.cpus} ${mem} --pfile ref --read-freq pca.acount \\
        --score pca.eigenvec.allele 2 \$A1 header-read no-mean-imputation variance-standardize \\
        --score-col-nums \$P1-\$P2 --out ref_pcs
    rm -f prune.prune.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ .*//')
    END_VERSIONS
    """

    stub:
    """
    touch "ref.pgen"
    touch "ref.pvar"
    touch "ref.psam"
    touch "ref.afreq"
    touch "pca.acount"
    touch "pca.eigenvec.allele"
    touch "pca.cols"
    touch "ref_pcs.sscore"
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stub: true
    END_VERSIONS
    """
}
