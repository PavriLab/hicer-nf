process HICUP_MAP_READS {

    tag "$meta.id"
    memory = { genomeSizeType == 'large' ? 80.GB * task.attempt : 20.GB * task.attempt }

    input:
    tuple val(meta), path(reads)
    path(index, stageAs: 'bowtie2Index')
    val(genomeSizeType)

    output:
    tuple val(meta), file("${meta.id}/${meta.id}_1_2.pair.sam"),    emit: alignments
    tuple val(meta), file("${meta.id}/*summary*.txt"),              emit: reports

    shell:
    println(index)
    '''
    mkdir !{meta.id}
    hicup_mapper \
        --outdir !{meta.id} \
        --threads !{task.cpus} \
        --format Sanger \
        --index !{index} \
        --bowtie2 $(which bowtie2) \
        !{reads}
    '''
}
