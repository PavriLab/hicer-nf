process SAM_TO_BAM {

    tag "$meta.id"

    publishDir path: "${params.outputDir}/${meta.id}/bam/",
               mode: 'copy',
               overwrite: 'true',
               pattern: "*.bam",
               saveAs: {
                   file ->
                        if (file.endsWith(".bam")) file(file).getName()
                }

    input:
    tuple val(meta), file(alignments)

    output:
    file("${meta.id}.hicup.bam")

    script:
    """
    samtools view -bh -@ ${task.cpus} ${alignments} > ${meta.id}.hicup.bam
    """
}
