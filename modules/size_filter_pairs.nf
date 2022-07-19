process SIZE_FILTER_PAIRS {

    tag "$meta.id"

    input:
    tuple val(meta), file(alignments)

    output:
    tuple val(meta), file("${meta.id}_1_2.filt.sam"),                           emit: alignments
    tuple val(meta), file("*summary*.txt"), file("*.ditag_size_distribution"),  emit: reports

    script:
    """
    filterBySize.py \
        -i ${alignments} \
        --minDistance ${params.minMapDistance} \
        -o ${meta.id}_1_2.filt.sam
    """
}
