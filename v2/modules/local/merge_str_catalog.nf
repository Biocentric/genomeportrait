process MERGE_STR_CATALOG {
    tag "str_catalog"
    label 'process_single'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] ? params.report_sif : params.report_docker }"

    input:
    path pathogenic          // Illumina ExpansionHunter catalog (repeat-expansion disease loci)
    path codis               // assets/codis_long_loci_hg38.json (long CODIS/ESS loci)

    output:
    path "str_catalog_merged.json", emit: catalog
    path "versions.yml",            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // ExpansionHunter accepts exactly one --variant-catalog, and the stock Illumina catalog
    // holds only pathogenic repeat-expansion loci -- no CODIS marker at all. Concatenate the
    // two arrays so a single EH pass genotypes both sets.
    // NB the python body must start at column 0: <<'PY' preserves leading whitespace.
    """
    python3 - <<'PY'
import json

def load(p):
    with open(p) as fh:
        d = json.load(fh)
    return d if isinstance(d, list) else [d]

path = load("${pathogenic}")
codis = load("${codis}")
seen = {e.get("LocusId") for e in path}
merged = path + [e for e in codis if e.get("LocusId") not in seen]
dupes = [e.get("LocusId") for e in codis if e.get("LocusId") in seen]
if dupes:
    print("already in the pathogenic catalog, not re-added: " + ", ".join(dupes))
with open("str_catalog_merged.json", "w") as fh:
    json.dump(merged, fh, indent=2)
print("merged catalog: %d pathogenic + %d CODIS = %d loci"
      % (len(path), len(merged) - len(path), len(merged)))
PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """
}
