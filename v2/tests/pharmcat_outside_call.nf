#!/usr/bin/env nextflow
/*
    Run the PHARMCAT module alone with a CYP2D6 outside call, to check the module passes it
    to PharmCAT and that PharmCAT reports it.

        nextflow run tests/pharmcat_outside_call.nf -c tests/pharmcat_outside_call.config \
            --cram NA12878.chr22.cram --fasta GRCh38.fasta --outside_call '*3/*68+*4'
*/
include { PHARMCAT } from '../modules/local/pharmcat'

params.cram         = null
params.fasta        = null
params.outside_call = '*3/*68+*4'

// CYP2D6 diplotype labels anywhere in a PharmCAT report JSON
def cyp2d6Labels(o) {
    def out = []
    if (o instanceof Map) {
        if (o.geneSymbol == 'CYP2D6' && o.sourceDiplotypes) out << o.sourceDiplotypes[0].label
        o.values().each { v -> out.addAll(cyp2d6Labels(v)) }
    } else if (o instanceof List) {
        o.each { v -> out.addAll(cyp2d6Labels(v)) }
    }
    return out
}

workflow {
    def outside = file("${workflow.workDir}/outside_call.tsv")
    outside.text = "CYP2D6\t${params.outside_call}\n"
    def cram = file(params.cram, checkIfExists: true)
    PHARMCAT (
        channel.of([ [ id: cram.simpleName ], cram, file("${params.cram}.crai", checkIfExists: true), outside ]),
        file(params.fasta, checkIfExists: true),
        file("${params.fasta}.fai", checkIfExists: true)
    )
    PHARMCAT.out.json.map { _meta, json ->
        def rep = new groovy.json.JsonSlurper().parse(json)
        def found = cyp2d6Labels(rep)
        def got = found ? found[0] : 'none'
        log.info "CYP2D6 in the PharmCAT report: ${got} (outside call: ${params.outside_call})"
        if (got != params.outside_call) error "PharmCAT did not report the outside call"
    }
}
