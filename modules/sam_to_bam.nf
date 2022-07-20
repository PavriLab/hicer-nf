process SAM_TO_BAM {

    tag "$meta.id"

    publishDir path: "${params.outputDir}/${meta.id}/bam/",
               mode: 'copy',
               overwrite: 'true',
               pattern: "*.bam",
               saveAs: {
                   outfile ->
                        if (outfile.endsWith(".bam")) file(outfile).getName()
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
