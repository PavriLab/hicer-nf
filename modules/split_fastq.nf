process SPLIT_FASTQ {

    tag "$meta.id"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}/${meta.id}_*"), emit: reads

    script:
    def numberOfLinesPerSplit = params.readsPerSplit.toInteger() * 4

    """
    mkdir ${meta.id}

    # ampersand at the end of the line lets linux execute the commands in parallel
    zcat ${reads[0]} | \
    split -l ${numberOfLinesPerSplit} -a 4 --additional-suffix _1.fq - ${meta.id}/${meta.id}_ &

    zcat ${reads[1]} | \
    split -l ${numberOfLinesPerSplit} -a 4 --additional-suffix _2.fq - ${meta.id}/${meta.id}_

    wait
    """
}
