process HICUP_DEDUPLICATE_PAIRS {

    tag "$meta.id"

    input:
    tuple val(meta), path(alignments)

    output:
    tuple val(meta), file("${meta.id}_1_2.dedup.sam"),  emit: alignments
    tuple val(meta), file("*summary*.txt"),             emit: reports

    shell:
    '''
    hicup_deduplicator \
        --outdir . \
        !{alignments}

    if [ ! !{params.re} ]
    then
        getCisFractions.py \
            -i !{meta.id}_1_2.dedup.sam \
            -r *summary*.txt
    fi
    '''
}
