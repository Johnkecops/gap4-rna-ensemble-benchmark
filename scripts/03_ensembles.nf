// Stage D: one process per (RNA, tool, beam). Run from repo root:
//   nextflow run scripts/03_ensembles.nf -with-trace results/nf_trace.tsv
// MERGE-RNA is fitted separately (scripts/merge_fit.py) because it fits jointly across RNAs.
nextflow.enable.dsl = 2
params.root  = "${projectDir}/.."
params.beams = [100, 200, 500, 0]

process FOLD {
    tag "${rid}:${tool}:${beam}"
    maxForks 4
    errorStrategy 'ignore'
    input:
    tuple val(rid), val(tool), val(beam)
    output:
    stdout
    script:
    """
    python ${params.root}/scripts/run_tool.py ${rid} ${tool} ${beam}
    """
}

workflow {
    ids = Channel.fromPath("${params.root}/data/processed/split.tsv")
                 .splitCsv(header: true, sep: '\t')
                 .map { it.rna_id }
    jobs = ids.flatMap { rid ->
        [[rid, 'vienna', -1]] +
        params.beams.collect { b -> [rid, 'linearpartition', b] } +
        params.beams.collect { b -> [rid, 'linearsampling', b] }
    }
    FOLD(jobs) | view
}
