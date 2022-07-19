process HICUP_DEDUPLICATE_PAIRS {

    tag "$meta.id"

    input:
    tuple val(meta), file(alignments)

    output:
    tuple val(meta) file("${resplitName}/${resplitName}_1_2.dedup.sam"), emit: alignments
    tuple val(resplitName), file("${resplitName}/*summary*.txt"),        emit: reports

    script:
    resplitName = alignments.getName() - ~/(_1_2\.filt\.sam)?$/

    """
    mkdir ${resplitName}
    hicup_deduplicator \
        --outdir ${resplitName} \
        ${alignments}

    if [ ! ${params.re} ]
    then
        getCisFractions.py \
            -i ${resplitName}/${resplitName}_1_2.dedup.sam \
            -r ${resplitName}/*summary*.txt
    fi
    """
}