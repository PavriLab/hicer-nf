process MAKE_BASE_MATRIX {

    tag "$meta.id"

    input:
    tuple val(meta), path(pairs)
    file(genomeSizes)
    val(genomeName)
    val(baseResolution)

    output:
    tuple val(meta), file("${meta.id}_base.cool"), emit: matrix

    script:
    """
    # making sure chromosomes are sorted semantically to comply with higlass
    sort -k1,1 -V ${genomeSizes} > chromSizes.sort.tsv

    cooler cload pairs \
        --assembly ${genomeName} \
        -c1 2 -p1 3 -c2 4 -p2 5 \
        chromSizes.sort.tsv:${baseResolution} \
        ${pairs[0]} \
        ${meta.id}_base.cool
    """
}
