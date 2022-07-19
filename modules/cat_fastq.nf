process CAT_FASTQ {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads, stageAs: "input*/*")

    output:
    tuple val(meta), path("*.merged.fastq.gz"), emit: reads

    script:
    def readList = reads.collect( it.toString() )
    def read1 = []
    def read2 = []
    readList.eachWithIndex{ v, ix -> ( ix & 1 ? read2 : read1 ) << v }

    """
    cat ${read1.join(' ')} > ${meta.id}_1.merged.fastq.gz
    cat ${read2.join(' ')} > ${meta.id}_2.merged.fastq.gz
    """
}
