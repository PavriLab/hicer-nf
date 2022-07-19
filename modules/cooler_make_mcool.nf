process MAKE_MCOOL {

    tag "$meta.id"

    input:
    tuple val(meta), file(basematrix)

    output:
    tuple val(meta), file("${name}.mcool"), emit: matrix

    script:
    """
    cooler zoomify \
        -p ${task.cpus} \
        -r ${resolutions} \
        -o ${meta.id}.mcool \
        ${basematrix}
    """
}
