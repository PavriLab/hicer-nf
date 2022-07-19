process CAT_SAM {

    tag "$meta.id"

    input:
    tuple val(meta), path(alignments)

    output:
    tuple val(meta), file("${name}_1_2.dedup.sam"), emit: alignments

    script:
    """
    samtools view -H ${alignments[0]} > ${meta.id}_1_2.dedup.sam
    cat ${alignments} | grep -v '^@' >> ${meta.id}_1_2.dedup.sam
    """
}
