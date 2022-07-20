process SAM_TO_BAM {

    tag "$meta.id"

    input:
    tuple val(meta), file(alignments)

    output:
    file("${meta.id}.hicup.bam")

    script:
    """
    samtools view -bh -@ ${task.cpus} ${alignments} > ${meta.id}.hicup.bam
    """
}
