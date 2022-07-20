process RESPLIT_FILTERED_PAIRS {

    tag "$meta.id"

    input:
    tuple val(meta), path(alignments)

    output:
    tuple val(meta), file("${meta.id}/${meta.id}*.sam"), emit: alignments

    script:
    """
    mkdir ${meta.id}
    resplitByChromosome.py \
        -i ${alignments} \
        -o ${meta.id}/${meta.id}
    """
}
