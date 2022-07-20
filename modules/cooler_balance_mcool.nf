process BALANCE_MCOOL {

    tag "$meta.id"

    input:
    tuple val(meta), file(mcool)

    output:
    tuple val(meta), file("${mcool}"), emit: matrix

    script:
    """
    balanceMultiCooler.py -m ${mcool} -p ${task.cpus}
    """
}
