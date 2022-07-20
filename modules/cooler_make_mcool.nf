process MAKE_MCOOL {

    tag "$meta.id"

    input:
    tuple val(meta), file(basematrix)
    val(resolutions)

    output:
    tuple val(meta), file("${meta.id}.mcool"), emit: matrix

    script:
    """
    cooler zoomify \
        -p ${task.cpus} \
        -r ${resolutions} \
        -o ${meta.id}.mcool \
        ${basematrix}
    """
}
