process BALANCE_MCOOL {

    tag "$meta.id"

    publishDir  path: "${params.outputDir}/${meta.id}/matrices/",
                mode: 'copy',
                overwrite: 'true',
                pattern: "${meta.id}.mcool"

    input:
    tuple val(meta), file(mcool)

    output:
    tuple val(meta), file("${mcool}"), emit: matrix

    script:
    """
    balanceMultiCooler.py -m ${mcool} -p ${task.cpus}
    """
}
