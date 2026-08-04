process YLEAF {
    tag "$meta.id"
    label 'process_medium'

    // The published yleaf biocontainer has no data/ dir and no samtools/bcftools, so it cannot
    // run at all; use the image the pipeline builds itself (BUILD_YLEAF_CONTAINER).
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? params.yleaf_sif : params.yleaf_docker }"

    input:
    tuple val(meta), path(cram), path(crai)
    path  fasta
    path  fai
    path  yleaf_data     // data/ dir missing from the conda package (DOWNLOAD_YLEAF_DATA)
    val   ready          // gate on BUILD_YLEAF_CONTAINER

    output:
    tuple val(meta), path("*.y_haplogroup.txt"), emit: hg
    path "${meta.id}_yleaf",                     emit: results, optional: true
    path "versions.yml",                         emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Subset to chrY first: Yleaf only looks at the Y, and running it straight off a
    // whole-genome CRAM makes samtools idxstats/mpileup crawl (~4 min for idxstats alone on a
    // 35 GB CRAM). Extracting chrY takes ~30 s for a ~0.3 GB BAM, after which Yleaf finishes
    // in seconds.
    // The conda yleaf package has no data/ dir, and its config.txt (which points Yleaf at a
    // reference genome) sits on the read-only container filesystem, so assemble a writable
    // copy of the package here: package .py files + the fetched data/ + a config.txt pointing
    // at the staged GRCh38 FASTA (avoiding a redundant 3 GB reference download).
    """
    YCHR=\$(samtools view -H -T ${fasta} ${cram} | awk -F'\\t' '/^@SQ/{for(i=1;i<=NF;i++){if(\$i=="SN:chrY"){print "chrY"} else if(\$i=="SN:Y"){print "Y"}}}' | head -1)
    [ -n "\$YCHR" ] || YCHR=chrY
    samtools view -@ ${task.cpus} -b -T ${fasta} -o ${meta.id}.bam ${cram} "\$YCHR"
    samtools index -@ ${task.cpus} ${meta.id}.bam

    # Remove any stale ./yleaf FIRST: python puts the cwd on sys.path, so a leftover copy
    # would shadow the installed package and resolve PKG to itself.
    rm -rf ./yleaf
    PKG=\$(python3 -c 'import yleaf, os; print(os.path.dirname(yleaf.__file__))')
    cp -r "\$PKG" ./yleaf
    rm -rf ./yleaf/data
    cp -r ${yleaf_data} ./yleaf/data

    cat > ./yleaf/config.txt <<CFG
full hg19 genome fasta location =
hg19 chromosome Y fasta location =
full hg38 genome fasta location = \$(readlink -f ${fasta})
hg38 chromosome Y fasta location =
CFG

    PYTHONPATH="\$PWD" python3 -m yleaf.Yleaf \\
        -bam ${meta.id}.bam -rg hg38 -o ${meta.id}_yleaf -t ${task.cpus} -force \\
        || echo "Yleaf exited non-zero" >&2

    # Collect the prediction. A bare `find ... || fallback` does NOT work here: find exits 0
    # when it matches nothing, so the fallback would never fire and the process would finish
    # with no output file at all (the original "Missing output file(s)" failure).
    HG=\$(find ${meta.id}_yleaf -name 'hg_prediction.hg' -o -name '*.hg' 2>/dev/null | head -1)
    if [ -n "\$HG" ] && [ -s "\$HG" ]; then
        cp "\$HG" ${meta.id}.y_haplogroup.txt
    else
        printf 'Sample_name\\tHg\\tNote\\n%s\\tNA\\tno_Y_haplogroup_called\\n' "${meta.id}" > ${meta.id}.y_haplogroup.txt
    fi
    rm -f ${meta.id}.bam ${meta.id}.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yleaf: \$(python3 -c 'import yleaf; print(getattr(yleaf, "__version__", "3.2.1"))' 2>/dev/null || echo 3.2.1)
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
