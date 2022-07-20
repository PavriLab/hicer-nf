process HICUP_DEDUPLICATE_PAIRS {

    tag "$meta.id"

    input:
    tuple val(meta), path(alignments)

    output:
    tuple val(meta), file("${resplitName}_1_2.dedup.sam"), emit: alignments
    tuple val(resplitName), file("*summary*.txt"),         emit: reports

    shell:
    def resplitName = alignments.name - ~/(_1_2\.filt\.sam)?$/

    '''
    hicup_deduplicator \
        --outdir . \
        !{alignments}

    if [ ! !{params.re} ]
    then
        getCisFractions.py \
            -i !{resplitName}_1_2.dedup.sam \
            -r *summary*.txt
    fi
    '''
}
